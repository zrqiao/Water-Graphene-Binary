//
// Created by utena on 17-8-23.
//
#include "Environment& Algorithm.hpp"
#include "amber_netcdf.hpp"
#include "amber_parm_1.0.hpp"
#include "vector_calc.h"
#define z_axis_modify  2
#define deviation_y_center 40
#define deviation_x_center 40
#define start_nc 0
#define end_nc  0
#define z_axis_num 200
#define C_NUM 3772
//~ #define name_nc "water_ion_graphene_10a5"
int main(int argc, char* argv[])
{
    CreatDir("z_distribution");
    char dist[16];
    char frcmod[16];
    sprintf(dist,"%s",argv[1]);
    sprintf(frcmod,"%s",argv[2]);
    int heat_temp=std::stoi(argv[3]);
    int equ_temp=std::stoi(argv[4]);
    typedef std::vector<double>::size_type index;
    char file_prefix[128];
    sprintf(file_prefix,"../dis%s_%s_WATBOX/MD_HEAT%dK_EQ%dK/Equilibrium/density_dis%s_WATBOX",dist,frcmod,heat_temp,equ_temp,dist);
    char name_parm7[128];
    sprintf(name_parm7, "%s.parm7", file_prefix);
    std::vector<int *> jump_coor;
    std::ofstream outfile;
    char name_output[128];
    sprintf(name_output, "z_distribution/%s.dat", dist);
    outfile.open(name_output);
    std::cout << "program to calculate the Z_axis distribution" << "\n" << std::endl;
    std::cout << "calculate XYZ_limit" << std::endl;
    std::vector<double> z_axis_count(z_axis_num, 0);
    double C_z_coor_sum_1 = 0, C_z_coor_sum_2 = 0;
    double C_z_coor_average_1, C_z_coor_average_2;
    double C_x_coor_sum = 0, C_y_coor_sum = 0;
    double C_x_coor_center, C_y_coor_center;
    double Z_UP, Z_DOWN, Y_UP, Y_DOWN, X_UP, X_DOWN;
    amber_parm parm_nam(name_parm7);
    char name_nc[128];
    sprintf(name_nc, "%s_NVT.nc", file_prefix);
    nctraj data_nc(name_nc);
    for (index C_index = 0; C_index != C_NUM; ++C_index) {
        C_z_coor_sum_1 += data_nc.atom_coordinate(0, C_index)[2];
        C_z_coor_sum_2 += data_nc.atom_coordinate(0, (C_NUM + C_index))[2];
        C_y_coor_sum += data_nc.atom_coordinate(0, C_index)[1];
        C_x_coor_sum += data_nc.atom_coordinate(0, C_index)[0];
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
    std::cout << "dz: " << dz << std::endl;

    amber_parm parm_name(name_parm7);
    sprintf(name_nc, "%s_NVT.nc", file_prefix);
    nctraj nc_data(name_nc);
    int total_frame = nc_data.frames_number();
    std::vector<index> O_WAT_id = parm_name.id_by_type("OW");
    int frame = 0;
    int WATcount=0;
    std::vector<double> o_coor;

    while (frame < total_frame) {
        for (index i = 0; i != O_WAT_id.size(); ++i) {
            o_coor = nc_data.atom_coordinate(frame, O_WAT_id[i]);
            if (o_coor[2] < Z_UP && o_coor[2] > Z_DOWN && o_coor[0] < X_UP && o_coor[0] > X_DOWN &&
                o_coor[1] < Y_UP && o_coor[1] > Y_DOWN) {
                z_axis_count[floor((o_coor[2] - Z_DOWN) / dz)] += 1;
                WATcount+=1;
            }
        }
        frame += 1;
    }
    for (int i = 0; i != z_axis_count.size(); ++i) {
        outfile << std::setw(15) << i*dz-(Z_UP-Z_DOWN)/2 << std::setw(15) << z_axis_count[i]/(WATcount)  << std::endl;
        std::cout << std::setw(15) << i*dz -(Z_UP-Z_DOWN)/2<< std::setw(15) << z_axis_count[i]/(WATcount)  << std::endl;
    }
    outfile.close();

    return 0;
}






