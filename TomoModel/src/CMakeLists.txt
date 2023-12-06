find_package(LAPACK REQUIRED)
list(APPEND LIBS ${LAPACK_LIBRARIES})

add_library(TomoModelModules
  nrtype.f90
  nrutil.f90
  module_util.f90
  module_function.f90
  module_fourier.f90
  module_spline.f90
  module_model.f90
  rspln.f
  rsple.f
  prott.f
  rotmx2.f
  lagrange_poly.f90
  gll_library.f90)

target_link_libraries(TomoModelModules PRIVATE ${LIBS})

add_executable(ReadS20RTS ReadS20RTS.f90)
target_link_libraries(ReadS20RTS TomoModelModules)

		 
