#include <vtkGenericDataObjectReader.h>
#include <vtkStructuredPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkDataArray.h>
#include <vtkPointSource.h>
#include <vtkOctreePointLocator.h>
#include <vtkDataSetCollection.h>
#include <vtkMath.h>
#include <string>
#include <iostream>
#include <sstream>

using namespace std;

enum DATATYPE
{
	SCALAR,
	VECTOR,
};

int main(int argc, char *argv[])
{
	// Ensure a filename was specified
	std::string inputFilename;
	if (argc > 2)
	{
		std::cerr << "Usage: " << argv[0] << " InputFilename" << endl;
		return EXIT_FAILURE;
	}
	else if (argc == 1)
		inputFilename = "out";
	else // Get the filename from the command line
		inputFilename = argv[1];

	int fileNum = 1;
	string inputDir = "InputData/";
	string ouputDir = "OutputData/";
	int dim = 64;

	for (int i = 0; i < fileNum; i++) {
		stringstream ss;
		ss << inputDir << inputFilename << "." << i + 1 << ".vtk" << endl;
		string fileName;
		ss >> fileName;
		// Get all data from the file
		vtkSmartPointer<vtkGenericDataObjectReader> reader =
			vtkSmartPointer<vtkGenericDataObjectReader>::New();
		reader->SetFileName(fileName.c_str());
		reader->SetReadAllScalars(true);
		reader->SetReadAllVectors(true);
		reader->Update();
		if (reader->IsFileUnstructuredGrid())
		{
			ss.clear();
			ss << ouputDir << inputFilename << "_" << dim << "_" << i << "_" << "structured" << ".vtk" << endl;
			std::string outputFilename;
			ss >> outputFilename;
			ofstream out(outputFilename);

			//write vtk head
			out << "# vtk DataFile Version 3.0" << endl;
			out << "unstructured grid data from " << inputDir << inputFilename << endl;
			out << "ASCII" << endl;
			out << "DATASET STRUCTURED_GRID" << endl;

			//output parameter
			out << "DIMENSIONS " << dim << " " << dim << " " << dim << endl;
			long pointsize = pow(dim, 3);
			out << "POINTS " << pointsize << " float" << endl;

			vtkUnstructuredGrid* output = reader->GetUnstructuredGridOutput();

			//output->GetAttributes();
			vtkFieldData* fd = output->GetAttributesAsFieldData(0);
			vtkPointData* pd = output->GetPointData();
			//Create the tree
			vtkSmartPointer<vtkOctreePointLocator> octree =
				vtkSmartPointer<vtkOctreePointLocator>::New();
			octree->SetDataSet(output);
			octree->BuildLocator();
			std::cout << "Number of points in tree: " << octree->GetDataSet()->GetNumberOfPoints() << std::endl;


			// out points positions
			double* bound = output->GetBounds();
			vector<double> datasetbox(bound, bound + 6);
			cout << "bounding box: ";
			for (auto b : datasetbox)
				cout << " " << b;
			cout << endl;

			for (int i = 0; i < dim; i++)
			{ //z
				double z = double((datasetbox[5] - datasetbox[4])*(i + 0.5) / dim + datasetbox[4]);
				for (int j = 0; j < dim; j++)
				{ //y
					double y = double((datasetbox[3] - datasetbox[2])*(j + 0.5) / dim + datasetbox[2]);
					for (int k = 0; k < dim; k++)
					{ //x
						double x = double((datasetbox[1] - datasetbox[0])*(k + 0.5) / dim + datasetbox[0]);
						out << x << " " << y << " " << z << endl;
					}
				}
			}

			//string vectorName = "vorticity";
			int pd_array_num = pd->GetNumberOfArrays();
			for (int aryId = 0; aryId < pd_array_num; aryId++)
			{
				//t cd_array_num = c

				const char* attributeName = pd->GetArrayName(aryId);

				if (attributeName == NULL)
					std::cout << "No such arry named: " << attributeName << std::endl;
				std::cout << "Array Name: " << attributeName << std::endl;

				if (!pd->HasArray(attributeName))
					std::cout << "No such arry named: " << attributeName << std::endl;
				vtkDataArray* dataflds = pd->GetArray(attributeName);

				DATATYPE type;
				float max = dataflds->GetMaxNorm();
				if (dataflds->GetNumberOfComponents() == 3)
					type = VECTOR;
				else if (dataflds->GetNumberOfComponents() == 1)
					type = SCALAR;
				else // ignore other data type
					continue;
				cout << "valid type: " << (type == SCALAR ? "sclar" : "vector") << endl;
				// output points attribute
				out << "POINT_DATA " << pointsize << endl;
				if (type == SCALAR)
					out << "SCALARS " << attributeName << " double" << endl;
				else
					out << "VECTORS " << attributeName << " double" << endl;

				for (int i = 0; i < dim; i++) { //z
					double z = double((datasetbox[5] - datasetbox[4])*(i + 0.5) / dim + datasetbox[4]);
					for (int j = 0; j < dim; j++) { //y
						double y = double((datasetbox[3] - datasetbox[2])*(j + 0.5) / dim + datasetbox[2]);
						for (int k = 0; k < dim; k++) //x
						{
							double x = double((datasetbox[1] - datasetbox[0])*(k + 0.5) / dim + datasetbox[0]);
							double point[3] = { x, y, z };
							vtkSmartPointer<vtkIdList> result =
								vtkSmartPointer<vtkIdList>::New();
							double radius = 2 * (datasetbox[1] - datasetbox[0]) / dim;
							// two interpolation methods
							// 1. interpolate by kNN points
							// 2. interpolate by points within radius

							//octree->FindClosestNPoints(8, point, result); //find 8 closest point
							octree->FindPointsWithinRadius(radius, point, result);

							//int closestId = octree->FindClosestPoint(point);
							float vectorres[3] = { 0, 0, 0 };
							float scalarres = 0;
							double accudist = 0;
							for (int n = 0; n < result->GetNumberOfIds(); n++) {
								double* referpt = output->GetPoint(result->GetId(n));
								double dis = vtkMath::Distance2BetweenPoints(referpt, point);
								double refersclar = 0;
								double* refervctor = nullptr;
								switch (type)
								{
								case SCALAR:
									refersclar = dataflds->GetTuple1(result->GetId(n));
									scalarres += refersclar * dis;
									break;
								case VECTOR:
									refervctor = dataflds->GetTuple3(result->GetId(n));
									for (int coord = 0; coord < 3; coord++)
										vectorres[coord] += refervctor[coord] * dis;
									break;
								default:
									break;
								}
								accudist += dis;
							}
							switch (type)
							{
							case SCALAR:
								if (accudist != 0)
									scalarres /= accudist;
								scalarres /= max;
								//if (scalarres > 0.01)
								out << scalarres << endl;
								break;
							case VECTOR:
								if (accudist != 0) {
									for (int coord = 0; coord < 3; coord++)
										vectorres[coord] /= accudist;
								}
								//vtkMath::Normalize(vectorres);
								out << vectorres[0] << " " << vectorres[1] << " " << vectorres[2] << endl;
								break;
							default:
								break;
							}
						}
					}
				}
			}
			out.close();
		}
	}
	return EXIT_SUCCESS;
}
