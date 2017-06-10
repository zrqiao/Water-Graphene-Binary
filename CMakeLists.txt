cmake_minimum_required(VERSION 3.6)
project(qzr)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_C_COMPILER g++)
set(SOURCE_FILES
        9a5/enthaply/amber_netcdf.hpp
        9a5/enthaply/amber_parm_1.0.hpp
        #9a5/cxx_diffusion.cpp
        9a5/Cavity.cpp
        #
        #18a5/z_distribution_wat.cpp
        #9a5/cavity_jump_correlation.cpp
        #9a5/find_jump_reverse.cpp
        #9a5/find_jump_reverse_start.cpp
        )
set(SOURCE_FILES1
        9a5/enthaply/amber_netcdf.hpp
        9a5/enthaply/amber_parm_1.0.hpp
        9a5/find_jump.cpp
        )
set(SOURCE_FILES2
        9a5/enthaply/amber_netcdf.hpp
        9a5/enthaply/amber_parm_1.0.hpp
        9a5/find_jump_reverse_start.cpp
        )

add_executable(qzr ${SOURCE_FILES})
add_executable(0a25ps ${SOURCE_FILES})
add_executable(qzr1 ${SOURCE_FILES1})
add_executable(qzr2 ${SOURCE_FILES2})

include_directories(/usr/local/include)
link_libraries(/usr/local/lib)
link_directories(/usr/local/lib)
target_link_libraries(qzr -lnetcdf_c++)
target_link_libraries(qzr1 -lnetcdf_c++)
target_link_libraries(qzr2 -lnetcdf_c++)
target_link_libraries(0a25ps -lnetcdf_c++)