//
// Created by utena on 17-7-10.
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
#define dt 1
#define frame_extend 0
#define stable_width_cutoff 0.24
#define dr 0.01
#define name_parm7 "density_dis9a5.parm7"
#define t_com 1000
typedef std::vector<double>::size_type index;
int judge_layer(std::vector<double> coor,double Z_DOWN,double Z_UP){
    if (coor[2] > (Z_DOWN + Z_UP) / 2) {
        return 1;
    }
    else {return -1;}
}
double calc_angle_z_axis(std::vector<double> ori){
    std::vector<double> z_axis={0,0,1};
    double angle;
    angle = vector_calc::vector_angle(z_axis,ori);
    return 180.0*acos(angle)/vector_calc::PI;
}
int main() {
    static double H1_avernum_upperlayer [z_points]={0};
    static double H2_avernum_upperlayer [z_points]={0};
    static double H1_avernum_lowerlayer [z_points]={0};
    static double H2_avernum_lowerlayer [z_points]={0};
    static double H1_avernum_upperlayer_t [t_com*2]={0};
    static double H2_avernum_upperlayer_t [t_com*2]={0};
    static double H1_avernum_lowerlayer_t [t_com*2]={0};
    static double H2_avernum_lowerlayer_t [t_com*2]={0};
    static double H1_correlationfunction_t [t_com*2]={0};
    static double H2_correlationfunction_t [t_com*2]={0};
    static double Nor_distribution[t_com*2][90]={0};
    static double H2_averz_t [t_com*2]={0};
    static double H1_averz_t [t_com*2]={0};
    static double O_averz_t [t_com*2]={0};
    static double H2_averangle_t [t_com*2]={0};
    static double H1_averangle_t [t_com*2]={0};
    static double Dipole_averangle_t [t_com*2]={0};
    static double Normal_averangle_t [t_com*2]={0};
    static std::vector<double> total_wat(z_points,0);
    static std::vector<double> total_wat_t(t_com*2,0);
    double dz, Z_UP, Z_DOWN, Y_UP, Y_DOWN, X_UP, X_DOWN;
    std::vector<double> OH1_vec(3,0),OH2_vec(3,0),Dip_vec(3,0),Nor_vec(3,0);
    std::ifstream infile;
    infile.open("H-Bond/H_Bond_O_index_transition_donor_up_1a2_narrow");
    std::ifstream infile0;
    infile0.open("H-Bond/H_Bond_O_index_transition_donor_up_zeropoint_1a2_narrow");
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
    std::vector<double> O1_coor, O2_coor, H1O_coor, H2O_coor,H1_coor, H2_coor;
    double angleH1,angleH2,angleDip,angleNor;
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
    int jump_count=2;
    double num[6];
    double zero_num[6];
    int Target_ID,nc,frame,HB_ID1,HB_ID2;
    infile >> num[0] >> num[1] >> num[2] >> num[3] >> num[4] >> num[5];
    int dtp;
    while (!infile0.fail()) {
        infile0 >> zero_num[0] >> zero_num[1] >> zero_num[2] >> zero_num[3] >> zero_num[4] >> zero_num[5];
        jump_count=zero_num[0];
        int nc0=zero_num[2],frame0=zero_num[3];
        std::cout<<jump_count<<std::endl;
        int temp_H1=zero_num[4],temp_H2=zero_num[5];//H1 lower H2 upper
        while (!infile.fail()) {
            if (jump_count==num[0]) {
                nc = num[2];
                frame = num[3];
                Target_ID = num[1];
                HB_ID1 = num[4];
                HB_ID2 = num[5];
                if (nc==nc0){
                    dtp=t_com+(frame-frame0)/dt;
                }
                else if (nc<nc0){
                    sprintf(name_nc, "nc/density_dis9a5_%d.nc", nc);
                    nctraj nc_data(name_nc);
                    dtp = t_com-(nc_data.frames_number()-frame+frame0)/dt;
                }
                else if (nc>nc0){
                    sprintf(name_nc, "nc/density_dis9a5_%d.nc", nc0);
                    nctraj nc_data(name_nc);
                    dtp = t_com+(nc_data.frames_number()-frame0+frame)/dt;
                }
                sprintf(name_nc, "nc/density_dis9a5_%d.nc", nc);
                nctraj nc_data(name_nc);
                O1_coor = nc_data.atom_coordinate(frame, Target_ID);
                H1_coor = nc_data.atom_coordinate(frame, Target_ID+temp_H1);
                H2_coor = nc_data.atom_coordinate(frame, Target_ID+temp_H2);
                OH1_vec = vector_calc::vector_minus(H1_coor,O1_coor);
                OH2_vec = vector_calc::vector_minus(H2_coor,O1_coor);
                Dip_vec = vector_calc::vector_minus(vector_calc::vector_divide(vector_calc::vector_add(H1_coor,H2_coor),2),O1_coor);
                Nor_vec = vector_calc::cross_product(OH2_vec,OH1_vec);
                angleH1=calc_angle_z_axis(OH1_vec);
                angleH2=calc_angle_z_axis(OH2_vec);
                angleDip=calc_angle_z_axis(Dip_vec);
                std::vector<double> z_axis={0,0,1};
                angleNor=180.0*acos((fabs(vector_calc::vector_angle(z_axis,Nor_vec))))/vector_calc::PI;
                //if (dtp==t_com) {std::cout<<angleNor<<std::endl;}
                total_wat[int(round((O1_coor[2] - Z_DOWN) / dz))]++;
                total_wat_t[dtp]++;
                H1_averz_t[dtp]+=H1_coor[2];
                H2_averz_t[dtp]+=H2_coor[2];
                O_averz_t[dtp]+=O1_coor[2];
                H1_averangle_t[dtp]+=angleH1;
                H2_averangle_t[dtp]+=angleH2;
                Dipole_averangle_t[dtp]+=angleDip;
                Normal_averangle_t[dtp]+=angleNor;
                Nor_distribution[dtp][int(floor(angleNor))]++;
                if (num[temp_H1+3] != -1) {
                    H1O_coor = nc_data.atom_coordinate(frame, num[temp_H1+3]);
                    if (judge_layer(H1O_coor, Z_DOWN, Z_UP) == 1) {
                        H1_avernum_upperlayer[int(round((O1_coor[2] - Z_DOWN) / dz))]++;
                        H1_avernum_upperlayer_t[dtp]++;
                    } else if (judge_layer(H1O_coor, Z_DOWN, Z_UP) == -1) {
                        H1_avernum_lowerlayer[int(round((O1_coor[2] - Z_DOWN) / dz))]++;
                        H1_avernum_lowerlayer_t[dtp]++;
                        H1_correlationfunction_t[dtp]+=1;
                    }
                }
                if (num[temp_H2+3] != -1) {
                    H2O_coor = nc_data.atom_coordinate(frame, num[temp_H2+3]);
                    if (judge_layer(H2O_coor, Z_DOWN, Z_UP) == 1) {
                        H2_avernum_upperlayer[int(round((O1_coor[2] - Z_DOWN) / dz))]++;
                        H2_avernum_upperlayer_t[dtp]++;
                        H2_correlationfunction_t[dtp]+=1;
                    } else if (judge_layer(H2O_coor, Z_DOWN, Z_UP) == -1) {
                        H2_avernum_lowerlayer[int(round((O1_coor[2] - Z_DOWN) / dz))]++;
                        H2_avernum_lowerlayer_t[dtp]++;
                    }
                }
            }
            else if(jump_count<num[0]){ break;}
            if (! infile.fail()) {
                infile >> num[0] >> num[1] >> num[2] >> num[3] >> num[4] >> num[5];
            }
            //std::cout<<Target_ID<<std::setw(10)<<frame<<std::endl;
        }
    }
    std::ofstream outfile;
    outfile.open("jump/cutoff_1a2/H1H2_z_averHbondnum_cutoff1a2_up_test");//0 z 1 H1u 2 H1l 3 H2u 4 H2l
    std::ofstream outfile1;
    outfile1.open("jump/cutoff_1a2/H1H2_t_averHbondnum_cutoff1a2_up_test");//0 t 1 H1u 2 H1l 3 H2u 4 H2l
    std::ofstream outfile2;
    outfile2.open("jump/cutoff_1a2/H1H2_t_averz_cutoff1a2_up_test");//0 z 1 H1 2 H2 3 O
    std::ofstream outfile3;
    outfile3.open("jump/cutoff_1a2/H1H2_t_averangle_zaxis_cutoff1a2_up_test");//0 t 1 H1 2 H2 3 D 4 N
    std::ofstream outfile4;
    outfile4.open("jump/cutoff_1a2/H1H2_t_Hbond_layerexchange_correlationfunction_test");
    for(int i =1; i !=z_points; ++i)
    {
        outfile<<Z_DOWN+i*dz<<std::setw(12)<<H1_avernum_upperlayer[i]/total_wat[i]<<std::setw(12)<<H1_avernum_lowerlayer[i]/total_wat[i]
               <<std::setw(12)<<H2_avernum_upperlayer[i]/total_wat[i]<<std::setw(12)<<H2_avernum_lowerlayer[i]/total_wat[i]<<std::endl;
        //std::cout <<total_wat[i]<<std::setw(2);
    }
    for(int i =1; i !=t_com*2; ++i)
    {
        outfile1<<4*(i*dt-t_com*dt)<<std::setw(12)<<H1_avernum_upperlayer_t[i]/total_wat_t[i]<<std::setw(12)<<H1_avernum_lowerlayer_t[i]/total_wat_t[i]
                <<std::setw(12)<<H2_avernum_upperlayer_t[i]/total_wat_t[i]<<std::setw(12)<<H2_avernum_lowerlayer_t[i]/total_wat_t[i]<<std::endl;
        outfile2<<4*(i*dt-t_com*dt)<<std::setw(12)<<H1_averz_t[i]/total_wat_t[i]<<std::setw(12)<<H2_averz_t[i]/total_wat_t[i]
                <<std::setw(12)<<O_averz_t[i]/total_wat_t[i]<<std::endl;
        outfile3<<4*(i*dt-t_com*dt)<<std::setw(12)<<H1_averangle_t[i]/total_wat_t[i]<<std::setw(12)<<H2_averangle_t[i]/total_wat_t[i]
                <<std::setw(12)<<Dipole_averangle_t[i]/total_wat_t[i]<<std::setw(12)<<Normal_averangle_t[i]/total_wat_t[i]<<std::endl;
        outfile4<<4*(i*dt-t_com*dt)<<std::setw(12)<<H1_correlationfunction_t[i]/total_wat_t[i]<<std::setw(12)<<H2_correlationfunction_t[i]/total_wat_t[i] <<std::endl;
        //std::cout <<total_wat[i]<<std::setw(2);

    }
    std::ofstream outfile5;
    outfile5.open("jump/cutoff_1a2/H1H2_t_Nor_angledistribution_test");
    for (int j=0;j<90;j++){
        for(int i =1; i !=t_com*2; ++i){
            outfile5<<Nor_distribution[i][j]/total_wat_t[i]<<std::setw(12);
        }
        outfile5<<std::endl;
    }
    outfile.close();
    outfile1.close();
    outfile2.close();
    outfile3.close();
    outfile4.close();
    outfile5.close();
    infile.close();
    infile0.close();
    return 0;

}