//
// Created by utena on 17-8-23.
//

#include "amber_netcdf.hpp"
#include "amber_parm_1.0.hpp"
#include "vector_calc.h"
#define z_axis_modify  2
#define deviation_y_center 50
#define deviation_x_center 50
#define start_nc 1
#define end_nc  5
#define z_axis_num 380
#define C_NUM 3772
//~ #define name_nc "water_ion_graphene_10a5"
int WAT_NUM[4]={1806,1700,2443,2550};
int main()
{
    typedef std::vector<double>::size_type index;
    for (int dens_ID=0;dens_ID<4;dens_ID++) {
        char name_parm7[64];
        sprintf(name_parm7, "WAT_%d/density_dis9a5_WAT%d.parm7", WAT_NUM[dens_ID], WAT_NUM[dens_ID]);
        std::vector<int *> jump_coor;
        std::ofstream outfile;
        std::ofstream outfile2;
        char name_output[64];
        sprintf(name_output, "WAT_%d/density_distribution_z_normalized", WAT_NUM[dens_ID]);
        outfile.open(name_output);
        std::cout << "program to calculate the Z_axis distribution" <<std::setw(5)<<WAT_NUM[dens_ID]<< "\n" << std::endl;
        std::cout << "calculate XYZ_limit" << std::endl;
        std::vector<double> z_axis_dis(z_axis_num, 0);
        double C_z_coor_sum_1 = 0, C_z_coor_sum_2 = 0;
        double C_z_coor_average_1, C_z_coor_average_2;
        double C_x_coor_sum = 0, C_y_coor_sum = 0;
        double C_x_coor_center, C_y_coor_center;
        double Z_UP, Z_DOWN, Y_UP, Y_DOWN, X_UP, X_DOWN;

        for (int nc = start_nc; nc <= end_nc; nc++) {
            amber_parm parm_nam(name_parm7);
            char name_nc[64];
            sprintf(name_nc, "WAT_%d/density_dis9a5_WAT%d_%d.nc", WAT_NUM[dens_ID],WAT_NUM[dens_ID],nc);
            nctraj data_nc(name_nc);
            for (index C_index = 0; C_index != C_NUM; ++C_index) {
                C_z_coor_sum_1 += data_nc.atom_coordinate(0, C_index)[2];
                C_z_coor_sum_2 += data_nc.atom_coordinate(0, (C_NUM + C_index))[2];
                C_y_coor_sum += data_nc.atom_coordinate(0, C_index)[1];
                C_x_coor_sum += data_nc.atom_coordinate(0, C_index)[0];
            }
        }
        C_z_coor_average_1 = C_z_coor_sum_1/(C_NUM*(end_nc-start_nc+1));
        C_z_coor_average_2 = C_z_coor_sum_2/(C_NUM*(end_nc-start_nc+1));
        C_y_coor_center = C_y_coor_sum/(C_NUM*(end_nc-start_nc+1));
        C_x_coor_center = C_x_coor_sum/(C_NUM*(end_nc-start_nc+1));

        Z_UP = C_z_coor_average_2;
        Z_DOWN = C_z_coor_average_1;
        Y_UP = C_y_coor_center + deviation_y_center;
        Y_DOWN = C_y_coor_center - deviation_y_center;
        X_UP = C_x_coor_center + deviation_x_center;
        X_DOWN = C_x_coor_center - deviation_x_center;
        double dz=(Z_UP-Z_DOWN)/z_axis_num;

        std::cout << "Z_UP: " << Z_UP << std::endl;
        std::cout << "Z_DOWN: " << Z_DOWN << std::endl;
        std::cout << "X_UP: " << X_UP << std::endl;
        std::cout << "X_DOWN: " << X_DOWN << std::endl;
        std::cout << "Y_UP: " << Y_UP << std::endl;
        std::cout << "Y_DOWN: " << Y_DOWN << std::endl;
        int frame_sum = 0;
        for (int nc = start_nc; nc != end_nc + 1; ++nc) {
            amber_parm parm_name(name_parm7);
            char name_nc[64];
            sprintf(name_nc, "WAT_%d/density_dis9a5_WAT%d_%d.nc", WAT_NUM[dens_ID],WAT_NUM[dens_ID],nc);
            nctraj nc_data(name_nc);
            int total_frame = nc_data.frames_number();
            std::vector<index> O_WAT_id = parm_name.id_by_type("OW");
            int frame = 0;
            std::vector<double> o_coor;
            double step;
            while (frame < total_frame) {
                for (index i = 0; i != O_WAT_id.size(); ++i) {
                    o_coor = nc_data.atom_coordinate(frame, O_WAT_id[i]);
                    if (o_coor[2] < Z_UP && o_coor[2] > Z_DOWN && o_coor[0] < X_UP && o_coor[0] > X_DOWN &&
                        o_coor[1] < Y_UP && o_coor[1] > Y_DOWN) {
                        step = (Z_UP - Z_DOWN) / z_axis_num;
                        z_axis_dis[floor((o_coor[2] - Z_DOWN) / step)] += 1;
                    }
                }
                frame += 10;
                frame_sum++;
            }
        }
        for (index i = 0; i != z_axis_dis.size(); ++i) {
            outfile << std::setw(15) << i*dz-(Z_UP-Z_DOWN)/2 << std::setw(15) << z_axis_dis[i]/(frame_sum*WAT_NUM[dens_ID])  << std::endl;
            std::cout << std::setw(15) << i*dz-(Z_UP-Z_DOWN)/2 << std::setw(15) << z_axis_dis[i]/(WAT_NUM[dens_ID]*frame_sum)  << std::endl;
        }
        outfile.close();
    }

    return 0;
}






