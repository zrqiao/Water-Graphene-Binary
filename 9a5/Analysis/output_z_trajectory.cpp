
#include "amber_netcdf.hpp"
#include "amber_parm_1.0.hpp"
#include "vector_calc.h"
//~ #include <string.h>
//~ #include <fstream>
//~ #include <iostream>
//~ #include <stdlib.h>
//~ #include <stdio.h>

#define z_axis_modify  2
#define deviation_y_center 26
#define deviation_x_center 26
#define start_nc 70
#define end_nc  80
#define select_boundary 3.55
#define hbond_cutoff_up 3.5
#define hbond_cutoff_down 2.0
#define hbond_cutoff_angle_up 180.0
#define hbond_cutoff_angle_down 135.0
#define bond_length_number 80
#define bond_angle_number 60
#define z_points    380
#define max_sampling 160000
#define dt 1
#define name_parm7 "density_dis9a5.parm7"
static double Z_UP, Z_DOWN, Y_UP, Y_DOWN, X_UP, X_DOWN;
//~ #define name_nc "water_ion_graphene_10a5"
typedef std::vector<double>::size_type index;
int main() {
    //std::ifstream infile1;
    std::vector<double> total_wat_z(z_points,0);
    std::vector<double> total_wat_z1(z_points,0);
    std::ifstream infile;
    std::ofstream outfile1,outfile2;
    std::cout << "program to calculate distribution during transition" << "\n" << std::endl;
    char name_nc[64];
    sprintf(name_nc, "nc/density_dis9a5_%d.nc", start_nc);
    amber_parm parm_name(name_parm7);
    nctraj nc_data(name_nc);
    std::vector<index> O_WAT_id = parm_name.id_by_type("OW");
    std::vector<index> O_WAT_IN_C_id, select_wat_id;
    std::vector<double> O_coor,O_coor_initial;
    double C_z_coor_sum_1 = 0, C_z_coor_sum_2 = 0;
    double C_z_coor_average_1, C_z_coor_average_2;
    double C_x_coor_sum = 0, C_y_coor_sum = 0;
    double C_x_coor_center, C_y_coor_center;
    double Z_UP,Z_DOWN,Y_UP,Y_DOWN,X_UP,X_DOWN;
    for (int nc=start_nc;nc<=end_nc;nc++) {
        amber_parm parm_nam(name_parm7);
        char name_nc[64];
        sprintf(name_nc, "nc/density_dis9a5_%d.nc",nc);
        nctraj data_nc(name_nc);
        for (index C_index = 0; C_index != 1400; ++C_index) {
            C_z_coor_sum_1 += data_nc.atom_coordinate(0, C_index)[2];
            C_z_coor_sum_2 += data_nc.atom_coordinate(0, (1400 + C_index))[2];
            C_y_coor_sum += data_nc.atom_coordinate(0, C_index)[1];
            C_x_coor_sum += data_nc.atom_coordinate(0, C_index)[0];
        }
    }
    C_z_coor_average_1 = C_z_coor_sum_1/(1400*(end_nc-start_nc+1));
    C_z_coor_average_2 = C_z_coor_sum_2/(1400*(end_nc-start_nc+1));
    C_y_coor_center = C_y_coor_sum/(1400*(end_nc-start_nc+1));
    C_x_coor_center = C_x_coor_sum/(1400*(end_nc-start_nc+1));

    Z_UP = C_z_coor_average_2;
    Z_DOWN = C_z_coor_average_1;
    Y_UP = C_y_coor_center + deviation_y_center;
    Y_DOWN = C_y_coor_center - deviation_y_center;
    X_UP = C_x_coor_center + deviation_x_center ;
    X_DOWN = C_x_coor_center - deviation_x_center ;

    std::cout << "Z_UP: " << Z_UP <<std::endl;
    std::cout << "Z_DOWN: " << Z_DOWN <<std::endl;
    std::cout << "X_UP: " << X_UP <<std::endl;
    std::cout << "X_DOWN: " << X_DOWN <<std::endl;
    std::cout << "Y_UP: " << Y_UP <<std::endl;
    std::cout << "Y_DOWN: " << Y_DOWN <<std::endl;
    /*infile.open("jump/cutoff_1a2/recrossing_index_start_finish_down_4");
    while (!infile.fail()) {
        double num[4];
        infile >> num[0] >> num[1]>>num[2]>>num[3];
        int frame_r_st = num[1]-750;
        int frame_r_ed = num[2]+750;
        int Target_ID = num[0];

        int nc = frame_r_st / 10000;
        sprintf(name_nc, "nc/density_dis9a5_%d.nc", nc);
        nctraj nc_data(name_nc);
        std::vector<double> O1_coor, O2_coor, H1_coor, H2_coor;
        if (frame_r_ed-frame_r_st>50) {
            nc=start_nc;
            int total_frame = (start_nc) * 10000;
            int totframe_nc=nc_data.frames_number();
            for (int frame_r = frame_r_st; frame_r < frame_r_ed; frame_r += dt) {
                int nc = start_nc;
                int total_frame = start_nc * 10000;
                if (frame_r==frame_r_st || frame_r>=(total_frame+totframe_nc)){
                    while (total_frame+totframe_nc <= frame_r) {
                        total_frame += totframe_nc;
                        nc++;
                        sprintf(name_nc, "nc/density_dis9a5_%d.nc", nc);
                        nctraj nc_data(name_nc);
                        totframe_nc=nc_data.frames_number();
                    }
                }
                sprintf(name_nc, "nc/density_dis9a5_%d.nc", nc);
                nctraj nc_data(name_nc);
                index frame = frame_r - total_frame ;
                O1_coor = nc_data.atom_coordinate(frame, Target_ID);
                std::cout << (frame_r-frame_r_st)*4<< std::setw(10) << O1_coor[2] << std::endl;
            }
            std::cout<<nc<<std::setw(10)<<Target_ID<<std::endl;
        }

        */
    std::ofstream outfile;
    outfile.open("RMS.test");
        //char name_nc[64];
        int nc=start_nc;
        int total_frame = (start_nc) * 10000;
        int totframe_nc=10000;
        int frame_r_st=start_nc*totframe_nc;
        for (int frame_r=start_nc*totframe_nc;frame_r<(end_nc+1)*totframe_nc;frame_r+=dt) {
            if (frame_r==frame_r_st || frame_r>=(total_frame+totframe_nc)){
                while (total_frame+totframe_nc <= frame_r) {
                    total_frame += totframe_nc;
                    nc++;
                    sprintf(name_nc, "nc/density_dis9a5_%d.nc", nc);
                    nctraj nc_data(name_nc);
                    totframe_nc=nc_data.frames_number();
                }
            }
            sprintf(name_nc, "nc/density_dis9a5_%d.nc", nc);
            nctraj nc_data(name_nc);

            int frame = frame_r - total_frame ;
            O_coor = nc_data.atom_coordinate(frame, O_WAT_id[0]);
            if ((frame_r-start_nc*totframe_nc)%40000==0){
                O_coor_initial=O_coor;
            }
            double dz=vector_calc::vector_minus(O_coor,O_coor_initial)[2];
            double dx=vector_calc::vector_minus(O_coor,O_coor_initial)[0];
            double dy=vector_calc::vector_minus(O_coor,O_coor_initial)[1];
            outfile<<frame_r*0.004<<std::setw(20)<<std::sqrt(pow(dx,2)+pow(dy,2))<<std::setw(20)<<std::sqrt((pow(dz,2)))<<std::endl;
            //total_wat_z[int(round((O1_coor[2] - Z_DOWN) / ((Z_UP - Z_DOWN) / z_points)))]++;
        }
        /*outfile1.open("jump/cutoff_1a2/z_distribution_with_transition_TSText");
        for (int i=0;i<z_points;i++){
            outfile1<<Z_DOWN+i*(Z_UP-Z_DOWN)/z_points<<std::setw(10)<< total_wat_z[i]<<std::endl;
        }
        outfile1.close();*/
    outfile.close();
    //infile1.close();
    return 0;
}
