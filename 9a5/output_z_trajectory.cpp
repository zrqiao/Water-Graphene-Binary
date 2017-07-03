
#include "amber_netcdf.hpp"
#include "amber_parm_1.0.hpp"
#include "vector_calc.h"
//~ #include <string.h>
//~ #include <fstream>
//~ #include <iostream>
//~ #include <stdlib.h>
//~ #include <stdio.h>


#define z_axis_modify  2
#define deviation_y_center 10
#define deviation_x_center 10
#define start_nc 37
#define end_nc  100
#define select_boundary 3.55
#define hbond_cutoff_up 3.5
#define hbond_cutoff_down 2.0
#define hbond_cutoff_angle_up 180.0
#define hbond_cutoff_angle_down 135.0
#define bond_length_number 80
#define bond_angle_number 60
#define z_points    160
#define max_sampling 160000
#define dt 4
#define name_parm7 "density_dis9a5.parm7"
static double Z_UP, Z_DOWN, Y_UP, Y_DOWN, X_UP, X_DOWN;
//~ #define name_nc "water_ion_graphene_10a5"
typedef std::vector<double>::size_type index;
int judge_H_Bond(std::vector<double> O1_coor,std::vector<double> H_coor,std::vector<double> O2_coor){
    double distance = sqrt(pow((O1_coor[0] -
                                O2_coor[0]), 2) +
                           pow((O1_coor[1] -
                                O2_coor[1]), 2) +
                           pow((O1_coor[2] -
                                O2_coor[2]), 2));
    if (hbond_cutoff_down < distance && distance < hbond_cutoff_up) {
        double angle;
        std::vector<double> coor1, coor2;
        coor1 = vector_calc::vector_minus(O1_coor, H_coor);
        coor2 = vector_calc::vector_minus(O2_coor, H_coor);
        angle = vector_calc::vector_angle(coor1, coor2);

        if (angle < cos(hbond_cutoff_angle_down * vector_calc::PI / 180.0)) {
            return 1;
        }
    }
    return 0;
}

int main() {
    //std::ifstream infile1;
    std::vector<double> total_wat_z(z_points,0);
    std::vector<double> total_wat_z1(z_points,0);
    double bond_distribution[z_points][z_points]={0};
    double bond_distribution_H[z_points][z_points]={0};
    std::ifstream infile;
    std::ofstream outfile1,outfile2;
    static std::vector<double> cavity_frames;
    std::cout << "program to calculate distribution during transition" << "\n" << std::endl;
    static std::vector<int> cavity_frame;

    index frame=0;
    index frame0=0;
    char name_nc[64];
    sprintf(name_nc, "density_dis9a5_%d.nc", start_nc);
    amber_parm parm_name(name_parm7);
    nctraj nc_data(name_nc);
    std::vector<index> O_WAT_id = parm_name.id_by_type("OW");
    std::vector<index> O_WAT_IN_C_id, select_wat_id;
    std::vector<double> O_coor,O_coor_initial;

    infile.open("transition_path_index_start_finish_down");
    while (!infile.fail()) {
        double num[4];
        infile >> num[0] >> num[1]>>num[2]>>num[3];
        int frame_r_st = num[1];
        int frame_r_ed = num[2];
        int Target_ID = num[0];

        int nc = frame_r_st / 10000;
        sprintf(name_nc, "density_dis9a5_%d.nc", nc);
        nctraj nc_data(name_nc);
        double C_z_coor_sum_1 = 0, C_z_coor_sum_2 = 0;
        double C_z_coor_average_1, C_z_coor_average_2;
        double C_x_coor_sum = 0, C_y_coor_sum = 0;
        double C_x_coor_center, C_y_coor_center;


        for (index C_index = 0; C_index != 1400; ++C_index) {
            C_z_coor_sum_1 += nc_data.atom_coordinate(frame, C_index)[2];
            C_z_coor_sum_2 += nc_data.atom_coordinate(frame, (1400 + C_index))[2];
            C_y_coor_sum += nc_data.atom_coordinate(frame, C_index)[1];
            C_x_coor_sum += nc_data.atom_coordinate(frame, C_index)[0];
        }
        C_z_coor_average_1 = C_z_coor_sum_1 / 1400;
        C_z_coor_average_2 = C_z_coor_sum_2 / 1400;
        C_y_coor_center = C_y_coor_sum / 1400;
        C_x_coor_center = C_x_coor_sum / 1400;

        Z_UP = C_z_coor_average_2;
        Z_DOWN = C_z_coor_average_1;
        Y_UP = C_y_coor_center + deviation_y_center;
        Y_DOWN = C_y_coor_center - deviation_y_center;
        X_UP = C_x_coor_center + deviation_x_center;
        X_DOWN = C_x_coor_center - deviation_x_center;
        std::vector<double> O1_coor, O2_coor, H1_coor, H2_coor;
        if (frame_r_ed-frame_r_st>2000) {
            for (int frame_r = frame_r_st; frame_r < frame_r_ed; frame_r += dt) {
                int nc = start_nc;
                int total_frame = start_nc * 10000;
                while (total_frame <= frame_r) {
                    sprintf(name_nc, "density_dis9a5_%d.nc", nc);
                    nctraj nc_data(name_nc);
                    total_frame += nc_data.frames_number();
                    nc++;
                }
                nc--;
                sprintf(name_nc, "density_dis9a5_%d.nc", nc);
                nctraj nc_data(name_nc);
                index frame = frame_r - total_frame + nc_data.frames_number();
                O1_coor = nc_data.atom_coordinate(frame, Target_ID);
                std::cout << frame_r - frame_r_st << std::setw(10) << O1_coor[2] << std::endl;
            }
        }
        /*for (int frame_r=frame_r_st;frame_r<frame_r_ed;frame_r+=dt) {
            int nc = start_nc;
            int total_frame = start_nc * 10000;
            while (total_frame <= frame_r) {
                sprintf(name_nc, "density_dis9a5_%d.nc", nc);
                nctraj nc_data(name_nc);
                total_frame += nc_data.frames_number();
                nc++;
            }
            nc--;
            sprintf(name_nc, "density_dis9a5_%d.nc", nc);
            nctraj nc_data(name_nc);
            index frame = frame_r - total_frame + nc_data.frames_number();
            select_wat_id.clear();
            O1_coor = nc_data.atom_coordinate(frame, Target_ID);
            if (num[3]==-1){
                total_wat_z[int(floor((O1_coor[2] - Z_DOWN) / ((Z_UP - Z_DOWN) / z_points)))]++;
            }
            else if(num[3]==1){
                total_wat_z1[int(floor((O1_coor[2] - Z_DOWN) / ((Z_UP - Z_DOWN) / z_points)))]++;
            }
        }
        std::cout<<Target_ID<<std::endl;
        outfile1.open("z_distribution_with_transition_down");
        outfile2.open("z_distribution_with_transition_up");
        for (int i=0;i<z_points;i++){
            outfile1<<Z_DOWN+i*(Z_UP-Z_DOWN)/z_points<<std::setw(10)<< total_wat_z[i]<<std::endl;
            outfile2<<Z_DOWN+i*(Z_UP-Z_DOWN)/z_points<<std::setw(10)<< total_wat_z1[i]<<std::endl;
        }
        outfile1.close();
        outfile2.close();*/
    }

    infile.close();

    //infile1.close();
    return 0;
}
