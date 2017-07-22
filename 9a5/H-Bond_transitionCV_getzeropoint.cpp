//
// Created by utena on 17-7-11.
//

#include "amber_netcdf.hpp"
#include "amber_parm_1.0.hpp"
#include "vector_calc.h"
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
#define z_points    320
#define max_sampling 160000
#define max_time 100
#define dt 4
#define frame_extend 0
#define stable_width_cutoff 0.24
#define dr 0.01
#define name_parm7 "density_dis9a5.parm7"
typedef std::vector<double>::size_type index;
int judge_layer(std::vector<double> coor,double Z_DOWN,double Z_UP){
    if (coor[2] > (Z_DOWN + Z_UP) / 2) {
        return 1;
    }
    else {return -1;}
}
int main() {
    static double H1_avernum_upperlayer [z_points]={0};
    static double H2_avernum_upperlayer [z_points]={0};
    static double H1_avernum_lowerlayer [z_points]={0};
    static double H2_avernum_lowerlayer [z_points]={0};
    static std::vector<double> total_wat(z_points,0);
    double dz, Z_UP, Z_DOWN, Y_UP, Y_DOWN, X_UP, X_DOWN;
    std::ifstream infile;
    static std::vector<double> cavity_frames;
    std::cout << "program to calculate zeropoint of jumping" << "\n" << std::endl;
    static std::vector<int> cavity_frame;
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
    std::vector<double> O1_coor, O2_coor, H1_coor, H2_coor;
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
    dz=(Z_UP-Z_DOWN)/z_points;
    std::cout << "Z_UP: " << Z_UP <<std::endl;
    std::cout << "Z_DOWN: " << Z_DOWN <<std::endl;
    std::cout << "X_UP: " << X_UP <<std::endl;
    std::cout << "X_DOWN: " << X_DOWN <<std::endl;
    std::cout << "Y_UP: " << Y_UP <<std::endl;
    std::cout << "Y_DOWN: " << Y_DOWN <<std::endl;
    infile.open("H-Bond/H_Bond_O_index_transition_donor_up_1a2_TSText");
    std::ifstream infile0;
    infile0.open("jump/cutoff_1a2/recrossing_index_start_finish_down_4");
    int jump_count;
    double num[6],num0[6];
    infile >> num[0] >> num[1] >> num[2] >> num[3] >> num[4] >> num[5];
    jump_count=num[0];
    int Target_ID,nc,frame,HB_ID1,HB_ID2;
    std:: ofstream outfile;
    outfile.open("H-Bond/H_Bond_O_index_transition_donor_up_zeropoint_1a2_TST");
    std::vector<int> zeropoint_stat(6,0);
    int frame0=0,framee=0;
    int jump_id=0;
    while (!infile0.fail()) {
        infile0 >> num0[0] >> num0[1]>>num0[2]>>num0[3];
        int frame_r_st = num0[1];
        int frame_r_ed = num0[2];
        int Target_ID = num0[0];
        int Direction = num0[3];
        if (jump_id==jump_count) {
            int temp_H1 = 0, temp_H2 = 0;//H1 lower H2 upper
            int O1_ID=-1,O2_ID=-1;
            int nc0=-1,frame0=-1;
            char name_nc[64];
            nc = start_nc;
            int total_frame = (start_nc) * 10000;
            int totframe_nc = nc_data.frames_number();
            while (total_frame + totframe_nc <= frame_r_st) {
                total_frame += totframe_nc;
                nc++;
                sprintf(name_nc, "nc/density_dis9a5_%d.nc", nc);
                nctraj nc_data(name_nc);
                totframe_nc = nc_data.frames_number();
            }
            int frame_st = frame_r_st - total_frame;
            int nc_st = nc;
            nc = start_nc;
            total_frame = (start_nc) * 10000;
            totframe_nc = nc_data.frames_number();
            while (total_frame + totframe_nc <= frame_r_ed) {
                total_frame += totframe_nc;
                nc++;
                sprintf(name_nc, "nc/density_dis9a5_%d.nc", nc);
                nctraj nc_data(name_nc);
                totframe_nc = nc_data.frames_number();
            }
            int frame_ed = frame_r_ed - total_frame;
            //std::cout << frame_st << std::setw(10) << frame_ed << std::endl;
            int nc_ed = nc;
            bool fixed = false,fixed_HB1=false,fixed_HB2=false;
            while (!infile.fail()) {
                if (jump_count != num[0]) {
                    break;
                } else {
                    nc = num[2];
                    frame = num[3];
                    Target_ID = num[1];
                    HB_ID1 = num[4];
                    HB_ID2 = num[5];
                }
                if (!fixed && ((frame >= frame_st && frame <= frame_ed)||(nc_st<nc_ed && frame>=frame_st)||(nc_st<nc_ed && frame<=frame_ed))) {//判断是否为参照点
                    sprintf(name_nc, "nc/density_dis9a5_%d.nc", nc);
                    nctraj nc_data(name_nc);
                    O1_coor = nc_data.atom_coordinate(frame, Target_ID);
                    if (O1_coor[2] > (Z_DOWN + Z_UP) / 2 - 0.06 && O1_coor[2] < (Z_DOWN + Z_UP) / 2 + 0.07) {
                        nc0 = nc;
                        frame0 = frame;
                        fixed = true;
                        /*if (HB_ID1 * HB_ID2 < 0) {
                            if (HB_ID1 == -1 &&
                                judge_layer(nc_data.atom_coordinate(frame, HB_ID2), Z_DOWN, Z_UP) == -1) {
                                temp_H1 = 2;
                                temp_H2 = 1;
                                fixed = true;
                            } else if (HB_ID2 == -1 &&
                                       judge_layer(nc_data.atom_coordinate(frame, HB_ID1), Z_DOWN, Z_UP) == -1) {
                                temp_H1 = 1;
                                temp_H2 = 2;
                                fixed = true;
                            } else if (HB_ID1 == -1 &&
                                       judge_layer(nc_data.atom_coordinate(frame, HB_ID2), Z_DOWN, Z_UP) == 1) {
                                temp_H1 = 1;
                                temp_H2 = 2;
                                fixed = true;
                            } else if (HB_ID2 == -1 &&
                                       judge_layer(nc_data.atom_coordinate(frame, HB_ID1), Z_DOWN, Z_UP) == 1) {
                                temp_H1 = 2;
                                temp_H2 = 1;
                                fixed = true;
                            }
                        } else if (HB_ID1 != -1 && HB_ID2 != -1) {
                            if (judge_layer(nc_data.atom_coordinate(frame, HB_ID1), Z_DOWN, Z_UP) == -1 &&
                                judge_layer(nc_data.atom_coordinate(frame, HB_ID2), Z_DOWN, Z_UP) == 1) {
                                temp_H1 = 1;
                                temp_H2 = 2;
                                fixed = true;
                            } else if (judge_layer(nc_data.atom_coordinate(frame, HB_ID2), Z_DOWN, Z_UP) == -1 &&
                                       judge_layer(nc_data.atom_coordinate(frame, HB_ID1), Z_DOWN, Z_UP) == 1) {
                                temp_H1 = 2;
                                temp_H2 = 1;
                                fixed = true;
                            }
                        }*/
                    }
                }
                if (!(fixed_HB1 && fixed_HB2) && ((frame >= frame_st && frame <= frame_ed)||(nc_st<nc_ed && frame>=frame_st)||(nc_st<nc_ed && frame<=frame_ed))) {//判断是否为参照点 HB1
                    sprintf(name_nc, "nc/density_dis9a5_%d.nc", nc);
                    nctraj nc_data(name_nc);
                    O1_coor = nc_data.atom_coordinate(frame, Target_ID);
                    if (O1_coor[2] > (Z_DOWN + Z_UP) / 2 - 0.2 && O1_coor[2] < (Z_DOWN + Z_UP) / 2 + 0.2) {
                        /*if (HB_ID1 * HB_ID2 < 0) {
                            if (HB_ID1 == -1 &&
                                judge_layer(nc_data.atom_coordinate(frame, HB_ID2), Z_DOWN, Z_UP) == -1) {
                                O1_ID=HB_ID2;
                                temp_H1 = 2;
                                temp_H2 = 1;
                                fixed_HB1 = true;
                            } else if (HB_ID2 == -1 &&
                                       judge_layer(nc_data.atom_coordinate(frame, HB_ID1), Z_DOWN, Z_UP) == -1) {
                                temp_H1 = 1;
                                temp_H2 = 2;
                                O1_ID=HB_ID1;
                                fixed_HB1 = true;
                            } else if (HB_ID1 == -1 &&
                                       judge_layer(nc_data.atom_coordinate(frame, HB_ID2), Z_DOWN, Z_UP) == 1) {
                                temp_H1 = 1;
                                temp_H2 = 2;
                                O2_ID=HB_ID2;
                                fixed_HB2 = true;
                            } else if (HB_ID2 == -1 &&
                                       judge_layer(nc_data.atom_coordinate(frame, HB_ID1), Z_DOWN, Z_UP) == 1) {
                                temp_H1 = 2;
                                temp_H2 = 1;
                                O2_ID=HB_ID1;
                                fixed_HB2 = true;
                            }
                        } else */if (HB_ID1 != -1 && HB_ID2 != -1) {
                            if (judge_layer(nc_data.atom_coordinate(frame, HB_ID1), Z_DOWN, Z_UP) == -1 &&
                                judge_layer(nc_data.atom_coordinate(frame, HB_ID2), Z_DOWN, Z_UP) == 1) {
                                temp_H1 = 1;
                                temp_H2 = 2;
                                O1_ID=HB_ID1;
                                O2_ID=HB_ID2;
                                fixed_HB1 = true;
                                fixed_HB2 = true;
                            } else if (judge_layer(nc_data.atom_coordinate(frame, HB_ID2), Z_DOWN, Z_UP) == -1 &&
                                       judge_layer(nc_data.atom_coordinate(frame, HB_ID1), Z_DOWN, Z_UP) == 1) {
                                temp_H1 = 2;
                                temp_H2 = 1;
                                O1_ID=HB_ID2;
                                O2_ID=HB_ID1;
                                fixed_HB1 = true;
                                fixed_HB2 = true;
                            }
                        }
                    }
                }
                if (!infile.fail()) {
                    infile >> num[0] >> num[1] >> num[2] >> num[3] >> num[4] >> num[5];
                }
                else { break; }
                //std::cout<<Target_ID<<std::setw(10)<<frame<<std::endl;
            }
            if (fixed && fixed_HB1 &&fixed_HB2) {
                outfile << jump_count << std::setw(10) << Target_ID << std::setw(10) << nc0 << std::setw(10)
                        << frame0 << std::setw(10) << temp_H1 << std::setw(10) << temp_H2
                        <<std::setw(10)<<O1_ID<<std::setw(10)<<O2_ID<< std::endl;
            }
            jump_count = num[0];
            std::cout << jump_id << std::endl;
        }
        jump_id++;
        //std::cout<<jump_id<<std::endl;
    }
    infile.close();
    outfile.close();
    return 0;

}