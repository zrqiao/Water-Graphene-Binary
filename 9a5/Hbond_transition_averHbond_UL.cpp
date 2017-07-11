//
// Created by utena on 17-7-11.
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
    double aver_Hbond_num_upperlayer_acceptor[z_points]={0};
    double aver_Hbond_num_lowerlayer_acceptor[z_points]={0};
    double aver_Hbond_num_upperlayer_donor[z_points]={0};
    double aver_Hbond_num_lowerlayer_donor[z_points]={0};
    static std::vector<double> total_wat_z(z_points,0);
    double dz, Z_UP, Z_DOWN, Y_UP, Y_DOWN, X_UP, X_DOWN;
    std::ifstream infile;
    static std::vector<double> cavity_frames;
    std::cout << "分别统计transitonpath上上下层acceptor donor的平均氢键数" << "\n" << std::endl;
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
    infile.open("H-Bond/H_Bond_O_index_during_transition_read");
    int count=0;
    int Target_ID=0,nc=0,frame=0,HB_ID=0;
    while (!infile.fail()) {
        double num[4];
        infile >> num[0] >> num[1]>>num[2]>>num[3];
        if (true) {
            if (!(num[0]==Target_ID&&nc==num[1]&&frame==num[2])){
                Target_ID = num[0];
                nc=num[1];
                frame = num[2];
                sprintf(name_nc, "nc/density_dis9a5_%d.nc", nc);
                nctraj nc_data(name_nc);
                O1_coor = nc_data.atom_coordinate(frame, Target_ID);
                total_wat_z[int(floor((O1_coor[2] - Z_DOWN) / dz))]++;
                count++;
            }
            HB_ID=num[3];
            std::cout<<nc<<' '<<frame<<' '<<Target_ID<<' '<<HB_ID<<std::endl;
            O2_coor = nc_data.atom_coordinate(frame, HB_ID);
            H1_coor = nc_data.atom_coordinate(frame, HB_ID + 1);
            H2_coor = nc_data.atom_coordinate(frame, HB_ID + 2);
            if (judge_layer(O2_coor, Z_DOWN, Z_UP) == 1) {
                aver_Hbond_num_upperlayer_acceptor[int(floor((O1_coor[2] - Z_DOWN) / ((Z_UP - Z_DOWN) / z_points)))]
                        += (judge_H_Bond(O1_coor, H1_coor, O2_coor) + judge_H_Bond(O1_coor, H2_coor, O2_coor));
            } else if (judge_layer(O2_coor, Z_DOWN, Z_UP) == -1) {
                aver_Hbond_num_lowerlayer_acceptor[int(floor((O1_coor[2] - Z_DOWN) / ((Z_UP - Z_DOWN) / z_points)))]
                        += (judge_H_Bond(O1_coor, H1_coor, O2_coor) + judge_H_Bond(O1_coor, H2_coor, O2_coor));
            }
            std::cout<<judge_H_Bond(O1_coor, H1_coor, O2_coor)<<judge_H_Bond(O1_coor, H2_coor, O2_coor);
            H1_coor = nc_data.atom_coordinate(frame, Target_ID + 1);
            H2_coor = nc_data.atom_coordinate(frame, Target_ID + 2);
            if (judge_layer(O2_coor, Z_DOWN, Z_UP) == 1) {
                aver_Hbond_num_upperlayer_donor[int(floor((O1_coor[2] - Z_DOWN) / ((Z_UP - Z_DOWN) / z_points)))]
                        += (judge_H_Bond(O1_coor, H1_coor, O2_coor) + judge_H_Bond(O1_coor, H2_coor, O2_coor));
            } else if (judge_layer(O2_coor, Z_DOWN, Z_UP) == -1) {
                aver_Hbond_num_lowerlayer_donor[int(floor((O1_coor[2] - Z_DOWN) / ((Z_UP - Z_DOWN) / z_points)))]
                        += (judge_H_Bond(O1_coor, H1_coor, O2_coor) + judge_H_Bond(O1_coor, H2_coor, O2_coor));
            }
            std::cout<<judge_H_Bond(O1_coor, H1_coor, O2_coor) <<judge_H_Bond(O1_coor, H2_coor, O2_coor)<<std::endl;
            //std::cout<<Target_ID<<std::setw(10)<<frame<<std::endl;
            if (count % 10000 == 0) {
                std::ofstream outfile;
                outfile.open("H-Bond/average_Hbond_num_ULlayer_transitionpath");//0 Z 1 AU 2 AL 3 DU 4 DL
                for (int i = 1; i != z_points; ++i) {
                    outfile << Z_DOWN + i * dz << std::setw(12)
                            << aver_Hbond_num_upperlayer_acceptor[i] / total_wat_z[i] << std::setw(12)
                            << aver_Hbond_num_lowerlayer_acceptor[i] / total_wat_z[i]
                            << std::setw(12) << aver_Hbond_num_upperlayer_donor[i] / total_wat_z[i] << std::setw(12)
                            << aver_Hbond_num_lowerlayer_donor[i] / total_wat_z[i] << std::endl;
                    //std::cout <<total_wat[i]<<std::setw(2);
                }
                outfile.close();
            }
        }
    }
    infile.close();
    return 0;

}