cmake_minimum_required(VERSION 3.0.2)
project(LADEL VERSION 0.1.0)

############################################################
################## Configure directories ###################
############################################################
get_directory_property(hasParent PARENT_DIRECTORY)
if(hasParent)
  # Assume the relevant directories have been set
else()
    SET(LIB_DIR ${CMAKE_CURRENT_BINARY_DIR}/../lib)

    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LIB_DIR})
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
endif()

############################################################
################### Compilation options ####################
############################################################
option (UNITTESTS "Perform unit testing" ON)
option (SIMPLE_COL_COUNTS "Simplify symbolic analysis (worse asymptotic time complexity)" OFF)

if (SIMPLE_COL_COUNTS)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DSIMPLE_COL_COUNTS")
endif (SIMPLE_COL_COUNTS) 

option (LONG "Use long integers for indexing" ON)
if (LONG)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DDLONG")
endif (LONG)

option (FLOAT "Use single precision instead of double" OFF)
if (FLOAT)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DDFLOAT")
endif (FLOAT)

if (UNIX)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -Wno-conversion -Wno-sign-conversion")
endif (UNIX)

option (AMD "Use approximate minimum degree ordering" ON)
if (AMD)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DDAMD")
endif (AMD)

option (PYTHON "Compile for python" OFF)

############################################################
#################### Build the library #####################
############################################################
SET(SRC
        src/ladel_global.c
        src/ladel_scale.c
        src/ladel_matvec.c
        src/ladel_matmat.c
        src/ladel_upper_diag.c
        src/ladel_etree.c
        src/ladel_postorder.c
        src/ladel_transpose.c
        src/ladel_col_counts.c
        src/ladel_ldl_symbolic.c
        src/ladel_permutation.c
        src/ladel_copy.c
        src/ladel_ldl_numeric.c
        src/ladel_pattern.c
        src/ladel.c
        src/ladel_rank1_mod.c
        src/ladel_row_mod.c
        src/ladel_debug_print.c
        src/ladel_add.c
        src/ladel_submatrix.c)


if (AMD)
    SET(AMD_SRC
            amd/Source/amd_1.c    
            amd/Source/amd_control.c   
            amd/Source/amd_global.c  
            amd/Source/amd_postorder.c   
            amd/Source/amd_valid.c
            amd/Source/amd_2.c    
            amd/Source/amd_defaults.c  
            amd/Source/amd_info.c    
            amd/Source/amd_post_tree.c
            amd/Source/amd_aat.c  
            amd/Source/amd_dump.c      
            amd/Source/amd_order.c   
            amd/Source/amd_preprocess.c
            amd/Source/SuiteSparse_config.c)
    SET(SRC ${SRC} ${AMD_SRC})
endif (AMD)

add_library(ladel ${SRC})

if (AMD)
    target_include_directories(ladel
        PUBLIC  ${CMAKE_CURRENT_LIST_DIR}/include 
                ${CMAKE_CURRENT_LIST_DIR}/amd/Include)
else ()
    target_include_directories(ladel
        PUBLIC  ${CMAKE_CURRENT_LIST_DIR}/include)

endif (AMD)

if (PYTHON)
    find_package(PythonInterp REQUIRED)
    # We don't really need the libs, but the headers.
    find_package(PythonLibs 3 REQUIRED)
    target_include_directories(ladel PUBLIC ${PYTHON_INCLUDE_DIRS})
    target_link_libraries(ladel ${PYTHON_LIBRARIES})
endif()
############################################################
#################### Build executables #####################
############################################################
if (UNITTESTS)
    add_subdirectory(test)
    include(CTest)
    enable_testing()
    add_test(NAME ladel_tester COMMAND $<TARGET_FILE:ladel_run_all_tests>)
endif()