cmake_minimum_required( VERSION 3.10 )

project( ex01 LANGUAGES CXX )

option ( ENABLE_LIKWID "ENABLE LIKWID MARKER"     ON )

find_package(BLAS)

add_executable(ex01_example ex01_example.cpp)


set_property(TARGET ex01_example PROPERTY CXX_STANDARD 17)

target_link_libraries( ex01_example PRIVATE BLAS::BLAS)

if ( ENABLE_LIKWID )
  find_library( LIKWID_LIB likwid PATH $ENV{LIKWID_LIBDIR})
  find_path( LIKWID_INCLUDE_DIR likwid.h PATH $ENV{LIKWID_INCDIR})
  if ( NOT LIKWID_LIB OR NOT LIKWID_INCLUDE_DIR )
    message(FATAL_ERROR "likwid not found. Help cmake with e.g. cmake .. -DCMAKE_PREFIX_PATH=/usr/local/likwid")
  endif()

  target_compile_definitions( ex01_example PRIVATE "-DUSE_LIKWID" )
  target_include_directories( ex01_example PRIVATE ${LIKWID_INCLUDE_DIR} )
  target_link_libraries(ex01_example PRIVATE ${LIKWID_LIB} )
endif()
