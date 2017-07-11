//
// Created by utena on 17-7-10.
//

#include "amber_netcdf.hpp"
#include "amber_parm_1.0.hpp"
#include "vector_calc.h"
#define deviation_y_center 10
#define deviation_x_center 10
#define start_nc 37
#define end_nc  52
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
    std::cout << "program to calculate average H bond number of H1 and H2" << "\n" << std::endl;
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
    infile.open("H-Bond/H_Bond_O_index_transition_donor_up");
    int jump_count=2;
    double num[6];
    infile >> num[0] >> num[1] >> num[2] >> num[3] >> num[4] >> num[5];
    int Target_ID,nc,frame,HB_ID1,HB_ID2;
    while (!infile.fail()) {
        int temp_H1=0,temp_H2=0;//H1 lower H2 upper
        bool fixed= false;
        while (!infile.fail()) {
            if (jump_count!=num[0]){break;}
            else {
                nc = num[2];
                frame = num[3];
                Target_ID = num[1];
                HB_ID1 = num[4];
                HB_ID2 = num[5];
            }
            sprintf(name_nc, "nc/density_dis9a5_%d.nc", nc);
            nctraj nc_data(name_nc);
            O1_coor = nc_data.atom_coordinate(frame, Target_ID);
            if (!fixed && O1_coor[2]>(Z_DOWN+Z_UP)/2-1.1){
                if (HB_ID1*HB_ID2<0){
                    if (HB_ID1==-1 && judge_layer(nc_data.atom_coordinate(frame,HB_ID2),Z_DOWN,Z_UP)==-1){
                        temp_H1=5;temp_H2=4;
                        fixed=true;
                    }
                    else if (HB_ID2==-1 && judge_layer(nc_data.atom_coordinate(frame,HB_ID1),Z_DOWN,Z_UP)==-1){
                        temp_H1=4;temp_H2=5;
                        fixed=true;
                    }
                }
                else if (HB_ID1!=-1 && HB_ID2!=-1 ){
                    if (judge_layer(nc_data.atom_coordinate(frame,HB_ID1),Z_DOWN,Z_UP)==-1 && judge_layer(nc_data.atom_coordinate(frame,HB_ID2),Z_DOWN,Z_UP)==1){
                        temp_H1=4;temp_H2=5;
                        fixed=true;
                    }
                    else if (judge_layer(nc_data.atom_coordinate(frame,HB_ID2),Z_DOWN,Z_UP)==-1 && judge_layer(nc_data.atom_coordinate(frame,HB_ID1),Z_DOWN,Z_UP)==1){
                        temp_H1=5;temp_H2=4;
                        fixed=true;
                    }
                }
            }
            if (fixed) {
                total_wat[int(round((O1_coor[2]-Z_DOWN)/dz))]++;
                if (num[temp_H1]!=-1){
                    H1_coor=nc_data.atom_coordinate(frame,num[temp_H1]);
                    if (judge_layer(H1_coor,Z_DOWN,Z_UP)==1){
                        H1_avernum_upperlayer[int(round((O1_coor[2]-Z_DOWN)/dz))]++;
                    }
                    else if (judge_layer(H1_coor,Z_DOWN,Z_UP)==-1){
                        H1_avernum_lowerlayer[int(round((O1_coor[2]-Z_DOWN)/dz))]++;
                    }
                }
                if (num[temp_H2]!=-1){
                    H2_coor=nc_data.atom_coordinate(frame,num[temp_H2]);
                    if (judge_layer(H2_coor,Z_DOWN,Z_UP)==1){
                        H2_avernum_upperlayer[int(round((O1_coor[2]-Z_DOWN)/dz))]++;
                    }
                    else if (judge_layer(H2_coor,Z_DOWN,Z_UP)==-1){
                        H2_avernum_lowerlayer[int(round((O1_coor[2]-Z_DOWN)/dz))]++;
                    }
                }
            }
            if (! infile.fail()) {
                infile >> num[0] >> num[1] >> num[2] >> num[3] >> num[4] >> num[5];
            }
            else {break;}
            //std::cout<<Target_ID<<std::setw(10)<<frame<<std::endl;
        }
        jump_count=num[0];
        std::cout<<jump_count<<std::endl;
        /*std::ofstream outfile1;
        outfile1.open("H-Bond/H1_z_averHbondnum_cutoff1a3_upperlayer");
        std::ofstream outfile2;
        outfile2.open("H-Bond/H1_z_averHbondnum_cutoff1a3_lowerlayer");
        for(int i =1; i !=z_points; ++i)
        {
            outfile1<<Z_DOWN+i*dz<<std::setw(12)<<H1_avernum_upperlayer[i]/total_wat[i]<<std::endl;
            //std::cout <<total_wat[i]<<std::setw(2);
        }
        //std::cout<<std::endl;
        for(int i =1; i !=z_points; ++i)
        {
            outfile2<<Z_DOWN+i*dz<<std::setw(14)<<H1_avernum_lowerlayer[i]/total_wat[i]<<std::endl;
        }
        outfile1.close();
        outfile2.close();
        std::ofstream outfile3;
        outfile3.open("H-Bond/H2_z_averHbondnum_cutoff1a3_upperlayer");
        std::ofstream outfile4;
        outfile4.open("H-Bond/H2_z_averHbondnum_cutoff1a3_lowerlayer");
        for(int i =1; i !=z_points; ++i)
        {
            outfile3<<Z_DOWN+i*dz<<std::setw(12)<<H2_avernum_upperlayer[i]/total_wat[i]<<std::endl;
        }
        for(int i =1; i !=z_points; ++i)
        {
            outfile4<<Z_DOWN+i*dz<<std::setw(14)<<H2_avernum_lowerlayer[i]/total_wat[i]<<std::endl;
        }
        outfile3.close();
        outfile4.close();*/
        std::ofstream outfile;
        outfile.open("H-Bond/H1H2_z_averHbondnum_cutoff1a3");//0 z 1 H1u 2 H1l 3 H2u 4 H2l
        for(int i =1; i !=z_points; ++i)
        {
            outfile<<Z_DOWN+i*dz<<std::setw(12)<<H1_avernum_upperlayer[i]/total_wat[i]<<std::setw(12)<<H1_avernum_lowerlayer[i]/total_wat[i]
                   <<std::setw(12)<<H2_avernum_upperlayer[i]/total_wat[i]<<std::setw(12)<<H2_avernum_lowerlayer[i]/total_wat[i]<<std::endl;
            //std::cout <<total_wat[i]<<std::setw(2);
        }
    }
    infile.close();
    return 0;

}