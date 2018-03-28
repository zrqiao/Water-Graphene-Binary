//
// Created by utena on 18-1-23.
//

#include "Environment& Algorithm.hpp"
#include "amber_netcdf.hpp"
#include "amber_parm_1.0.hpp"
#include "vector_calc.h"

#define D_VDW 3.800
#define wat_mw  18.01524
#define z_axis_modify  2
#define deviation_y_center 40
#define deviation_x_center 40
#define start_nc 0
#define end_nc  0
#define z_axis_num 200
#define C_NUM 3772

//~ #define name_nc "water_ion_graphene_10a5"
int main(int argc, char* argv[]) {
    char* dist_array[7]={"8a0","8a3","8a5","8a7","9a0","9a5","10a0"};
    typedef std::vector<double>::size_type index;
    std::ofstream outfile;
    char name_output[128];
    sprintf(name_output, "density_by_distance.dat");
    outfile.open(name_output);
    outfile<<"Distance"<<std::setw(15)<<"Density"<<std::endl;
    std::cout << "program to calculate the overall density" << "\n" << std::endl;
    char dist[16];
    char frcmod[16];

    for (int j=0; j<7; j++) {

        sprintf(dist, "%s", dist_array[j]);
        sprintf(frcmod, "%s", argv[2]);
        int heat_temp = std::stoi(argv[3]);
        int equ_temp = std::stoi(argv[4]);

        char file_prefix[128];
        sprintf(file_prefix, "../dis%s_%s_WATBOX/MD_HEAT%dK_EQ%dK/Equilibrium/density_dis%s_WATBOX", dist, frcmod,
                heat_temp, equ_temp, dist);
        char name_parm7[128];
        sprintf(name_parm7, "%s.parm7", file_prefix);
        amber_parm parm_nam(name_parm7);
        char name_nc[128];
        sprintf(name_nc, "%s_NVT.nc", file_prefix);
        nctraj data_nc(name_nc);

        std::cout << "calculate XYZ_limit" << std::endl;
        std::vector<double> z_axis_count(z_axis_num, 0);
        double C_z_coor_sum_1 = 0, C_z_coor_sum_2 = 0;
        double C_z_coor_average_1, C_z_coor_average_2;
        double C_x_coor_sum = 0, C_y_coor_sum = 0;
        double C_x_coor_center, C_y_coor_center;
        double Z_UP, Z_DOWN, Y_UP, Y_DOWN, X_UP, X_DOWN;

        for (index C_index = 0; C_index != C_NUM; ++C_index) {
            C_z_coor_sum_1 += data_nc.atom_coordinate(0, C_index)[2];
            C_z_coor_sum_2 += data_nc.atom_coordinate(0, (C_NUM + C_index))[2];
            C_y_coor_sum += data_nc.atom_coordinate(0, C_index)[1];
            C_x_coor_sum += data_nc.atom_coordinate(0, C_index)[0];
        }
        C_z_coor_average_1 = C_z_coor_sum_1 / (C_NUM * (end_nc - start_nc + 1));
        C_z_coor_average_2 = C_z_coor_sum_2 / (C_NUM * (end_nc - start_nc + 1));
        C_y_coor_center = C_y_coor_sum / (C_NUM * (end_nc - start_nc + 1));
        C_x_coor_center = C_x_coor_sum / (C_NUM * (end_nc - start_nc + 1));

        Z_UP = C_z_coor_average_2;
        Z_DOWN = C_z_coor_average_1;
        Y_UP = C_y_coor_center + deviation_y_center;
        Y_DOWN = C_y_coor_center - deviation_y_center;
        X_UP = C_x_coor_center + deviation_x_center;
        X_DOWN = C_x_coor_center - deviation_x_center;

        std::cout << "Z_UP: " << Z_UP << std::endl;
        std::cout << "Z_DOWN: " << Z_DOWN << std::endl;
        std::cout << "X_UP: " << X_UP << std::endl;
        std::cout << "X_DOWN: " << X_DOWN << std::endl;
        std::cout << "Y_UP: " << Y_UP << std::endl;
        std::cout << "Y_DOWN: " << Y_DOWN << std::endl;

        amber_parm parm_name(name_parm7);
        sprintf(name_nc, "%s_NVT.nc", file_prefix);
        nctraj nc_data(name_nc);
        int total_frame = nc_data.frames_number();
        std::vector <index> O_WAT_id = parm_name.id_by_type("OW");
        int frame = 0;
        int WATcount = 0;
        int frame_num = 0;
        std::vector<double> o_coor;

        while (frame < total_frame) {
            for (index i = 0; i != O_WAT_id.size(); ++i) {
                o_coor = nc_data.atom_coordinate(frame, O_WAT_id[i]);
                if (o_coor[2] < Z_UP && o_coor[2] > Z_DOWN && o_coor[0] < X_UP && o_coor[0] > X_DOWN &&
                    o_coor[1] < Y_UP && o_coor[1] > Y_DOWN) {
                    WATcount += 1;
                }
            }
            frame += 100;
            frame_num++;
        }

        double density=(WATcount/frame_num)/(4*deviation_x_center*deviation_y_center*(Z_UP-Z_DOWN-D_VDW))*wat_mw/0.622140857;
        outfile << (Z_UP - Z_DOWN) << std::setw(15) << density << std::endl;
        std::cout << (Z_UP - Z_DOWN) << std::setw(15) << density << std::endl;
    }

    outfile.close();

    return 0;
}






