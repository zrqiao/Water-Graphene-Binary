cmake_minimum_required(VERSION 3.6)
project(qzr)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_C_COMPILER g++)
set(SOURCE_FILES
        9a5/enthaply/amber_netcdf.hpp
        9a5/enthaply/amber_parm_1.0.hpp
        #9a5/cxx_diffusion.cpp
        #9a5/Cavity.cpp
        #9a5/find_jump.cpp
        #9a5/z_distribution_wat.cpp
        #9a5/cavity_jump_correlation.cpp
        #9a5/find_jump_reverse.cpp
        #9a5/find_jump_reverse_start.cpp
        #9a5/jump_time_distribution.cpp
        #9a5/jumptime_with_transition_test.cpp

        9a5/stabletime_distribution_during_transition
        )
add_executable(stable_time ${SOURCE_FILES})
target_link_libraries(stable_time -lnetcdf_c++)
add_executable(get_target_trajectory   9a5/enthaply/amber_netcdf.hpp
        9a5/enthaply/amber_parm_1.0.hpp 9a5/output_z_trajectory.cpp)
target_link_libraries(get_target_trajectory -lnetcdf_c++)
add_executable(format_translate   9a5/enthaply/amber_netcdf.hpp
        9a5/enthaply/amber_parm_1.0.hpp 9a5/format_translation_frame-nc.cpp)
target_link_libraries(format_translate -lnetcdf_c++)
set(SOURCE_Correlation
        9a5/enthaply/amber_netcdf.hpp
        9a5/enthaply/amber_parm_1.0.hpp
        #9a5/cxx_diffusion.cpp
        #9a5/find_jump.cpp
        #9a5/z_distribution_wat.cpp
        #9a5/cavity_jump_correlation.cpp
        #9a5/find_jump_reverse.cpp
        #9a5/find_jump_reverse_start.cpp
        #9a5/jump_time_distribution.cpp
        9a5/Hbond_x_z_transiton_correlation.cpp
        9a5/cavitysize_over_wat_with_transition.cpp)
add_executable(J-C_correlation ${SOURCE_Correlation})
target_link_libraries(J-C_correlation -lnetcdf_c++)
set(SOURCE_FILES1
        9a5/enthaply/amber_netcdf.hpp
        9a5/enthaply/amber_parm_1.0.hpp
        #

        9a5/Cavity.cpp
        #9a5/cavity_size_change-distribution.cpp
        #9a5/H-Bond_with_transition.cpp
        #water/Find_cavity_in_bulk_water.cpp 9a5/output_z_trajectory.cpp

        9a5/format_translation_frame-nc.cpp)
set(SOURCE_FILES2
        9a5/enthaply/amber_netcdf.hpp
        9a5/enthaply/amber_parm_1.0.hpp
        #9a5/find_jump_reverse_start.cpp
        #9a5/H_BOND_distribution.cpp
        #9a5/H-Bond_with_transition.cpp
        #9a5/H_BOND_distribution_O2.cpp
        9a5/H-Bond_by_z.cpp
        )


add_executable(jump ${SOURCE_FILES})
add_executable(cavity ${SOURCE_FILES1})
add_executable(H-Bond ${SOURCE_FILES2} )
add_executable(Bulk_water water/Find_cavity_in_bulk_water.cpp)

include_directories(/usr/local/include)
find_package(MPI)
include_directories(${MPI_INCLUDE_PATH})
link_libraries(/usr/local/lib)
link_directories(/usr/local/lib)

target_link_libraries(cavity -lnetcdf_c++ )
target_link_libraries(H-Bond -lnetcdf_c++)
target_link_libraries(jump -lnetcdf_c++)
target_link_libraries(Bulk_water -lnetcdf_c++)