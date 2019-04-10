
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

project(UnstructuredGridReader)

find_package(VTK COMPONENTS 
  vtkCommonCore
  vtkFiltersGeometry
  vtkIOLegacy
  vtkInteractionStyle
  vtkRenderingCore
  vtkRenderingFreeType
  vtkRenderingOpenGL2 QUIET)
if (NOT VTK_FOUND)
  message("Skipping UnstructuredGridReader: ${VTK_NOT_FOUND_MESSAGE}")
  return ()
endif()
message (STATUS "VTK_VERSION: ${VTK_VERSION}")
if (VTK_VERSION VERSION_LESS "8.90.0")
  # old system
  include(${VTK_USE_FILE})
  add_executable(UnstructuredGridReader MACOSX_BUNDLE UnstructuredGridReader.cxx )
  target_link_libraries(UnstructuredGridReader PRIVATE ${VTK_LIBRARIES})
else ()
  # include all components
  add_executable(UnstructuredGridReader MACOSX_BUNDLE UnstructuredGridReader.cxx )
  target_link_libraries(UnstructuredGridReader PRIVATE ${VTK_LIBRARIES})
  # vtk_module_autoinit is needed
  vtk_module_autoinit(
    TARGETS UnstructuredGridReader
    MODULES ${VTK_LIBRARIES}
    )
endif () 
