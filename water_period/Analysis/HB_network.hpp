//
// Created by utena on 17-10-15.
//

#include "Identification_of_grid_water.hpp"

using namespace std;

typedef std::vector<double>::size_type index0;
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
int Calc_HBnetwork(char* dist,int temperature,int start_nc,int end_nc) {
    CreatDir("Analysis/HB_Network");
    std::vector<double> OH1_vec(3, 0), OH2_vec(3, 0), Dip_vec(3, 0), Nor_vec(3, 0), O1H1_vec(3, 0), O2H2_vec(3, 0);
    char name_nc[64];
    std::vector <index0> O_WAT_IN_C_id, select_wat_id;
    double C_x_coor_sum = 0, C_y_coor_sum = 0;
    double C_x_coor_center, C_y_coor_center;
    std::vector<double> O_coor, O1_coor, O2_coor, H1O_coor, H2O_coor, H1_coor, H2_coor;
    double angleH1, angleH2, angleO1H1, angleO2H2, angleDip, angleNor;
    for (int nc = start_nc; nc <= end_nc; nc++) {
        char out_name[64];
        sprintf(out_name, "Analysis/HB_Network/nc_%02d.dat",nc);
        std::ofstream outfile;
        outfile.open(out_name);
        amber_parm parm_name(name_parm7);
        std::vector <index0> O_WAT_id = parm_name.id_by_type("OW");

        char name_nc[64];
        sprintf(name_nc, "density_dis%s_WAT1150_NVT_%d_%d.nc",dist,temperature,nc);
        nctraj nc_data(name_nc);
        int tot_frame=nc_data.frames_number();
        for (index0 frame = 0; frame < tot_frame; frame += 100) {
            for (index0 O_index = 0; O_index != O_WAT_id.size(); ++O_index) {
                C_y_coor_sum += nc_data.atom_coordinate(frame, O_WAT_id[O_index])[1];
                C_x_coor_sum += nc_data.atom_coordinate(frame, O_WAT_id[O_index])[0];
            }
        }
        C_y_coor_center = C_y_coor_sum / ( O_WAT_id.size()* (tot_frame/100));
        C_x_coor_center = C_x_coor_sum / ( O_WAT_id.size()* (tot_frame/100));

        double Y_UP = C_y_coor_center + y_length/2;
        double Y_DOWN = C_y_coor_center - y_length/2;
        double X_UP = C_x_coor_center + x_length/2;
        double X_DOWN = C_x_coor_center - x_length/2;
        X_DOWN=-0.4;
        Y_DOWN=-0.4;

        std::cout << "X_UP: " << X_UP << std::endl;
        std::cout << "X_DOWN: " << X_DOWN << std::endl;
        std::cout << "Y_UP: " << Y_UP << std::endl;
        std::cout << "Y_DOWN: " << Y_DOWN << std::endl;
        int Target_ID, frame, HB_ID1, HB_ID2;
        int x_points=ceil(x_length/select_boundary)+1;
        int y_points=ceil(y_length/select_boundary)+1;
        for (index0 frame = 0; frame < nc_data.frames_number(); frame += dt) {
            std::cout<<nc<<' '<<frame<<std::endl;
            std::vector < int > *xypoints[x_points+1];
            for (int i=0;i<x_points+1;i++){
                xypoints[i]=new std::vector<int>[y_points+1];
            }
            int check_connection[O_WAT_id.size()][O_WAT_id.size()]={0};
            int connection_data[O_WAT_id.size()][2];
            for (int i=0;i<O_WAT_id.size();i++){
                connection_data[i][0]=-1;
                connection_data[i][1]=-1;
            }
            if (frame % 100 == 0) std::cout << frame << std::endl;

            for (int i = 0; i != O_WAT_id.size(); i++) {//打格点
                //std::cout<<i<<std::endl;
                O_coor = nc_data.atom_coordinate(frame, O_WAT_id[i]);
                int x_location = floor((O_coor[0]-X_DOWN) / select_boundary);
                int y_location = floor((O_coor[1]-Y_DOWN) / select_boundary);
                xypoints[x_location][y_location].push_back(i);
                xypoints[x_location + 1][y_location].push_back(i);
                xypoints[x_location][y_location + 1].push_back(i);
                xypoints[x_location + 1][y_location + 1].push_back(i);
            }
            for (int i = 0; i < x_points; i++) {
                xypoints[i][0].insert(xypoints[i][0].end(), xypoints[i][y_points-1].begin(), xypoints[i][y_points-1].end());
            }
            for (int j = 0; j < y_points; j++) {
                xypoints[0][j].insert(xypoints[0][j].end(), xypoints[x_points-1][j].begin(), xypoints[x_points-1][j].end());
            }
            xypoints[0][0].insert(xypoints[0][0].end(), xypoints[x_points-1][y_points-1].begin(), xypoints[x_points-1][y_points-1].end());
            for (int i = 0; i < x_points; i++) {
                for (int j = 0; j < y_points; j++) {
                    for (vector<int>::iterator iter1 = xypoints[i][j].begin(); iter1 != xypoints[i][j].end(); iter1++) {
                        int wat1_i = *iter1;
                        for (vector<int>::iterator iter2 = xypoints[i][j].begin();
                             iter2 != xypoints[i][j].end(); iter2++) {
                            int wat2_i = *iter2;
                            if (check_connection[wat1_i][wat2_i]==0) {//避免重复计算
                                check_connection[wat1_i][wat2_i]=1;
                                O1_coor = nc_data.atom_coordinate(frame, O_WAT_id[wat1_i]);
                                O2_coor = nc_data.atom_coordinate(frame, O_WAT_id[wat2_i]);
                                H1_coor = nc_data.atom_coordinate(frame, O_WAT_id[wat1_i] + 1);
                                H2_coor = nc_data.atom_coordinate(frame, O_WAT_id[wat1_i] + 2);
                                if (judge_H_Bond(O1_coor, H1_coor, O2_coor)) {
                                    connection_data[wat1_i][0]=O_WAT_id[wat2_i];
                                }
                                if (judge_H_Bond(O1_coor, H2_coor, O2_coor)) {
                                    connection_data[wat1_i][1]=O_WAT_id[wat2_i];
                                }
                                std::vector<double> O_coor_period(3,0);
                                if (i==0){

                                    if (O2_coor[0]-O1_coor[0]>10){
                                        O_coor_period=O2_coor;
                                        O_coor_period[0]-=x_length;
                                    }
                                    else if (O2_coor[0]-O1_coor[0]<-10){
                                        O_coor_period=O2_coor;
                                        O_coor_period[0]+=x_length;
                                    }
                                    if (judge_H_Bond(O1_coor, H1_coor, O_coor_period)) {
                                        connection_data[wat1_i][0]=O_WAT_id[wat2_i];
                                        //std::cout<<O_WAT_id[wat1_i]<<"   "<<O_WAT_id[wat2_i]<<"  ";
                                    }
                                    if (judge_H_Bond(O1_coor, H2_coor, O_coor_period)) {
                                        connection_data[wat1_i][1]=O_WAT_id[wat2_i];
                                        //std::cout<<O_WAT_id[wat1_i]<<"   "<<O_WAT_id[wat2_i]<<"  ";
                                    }
                                }
                                if (j==0){
                                    if (O2_coor[1]-O1_coor[1]>10){
                                        O_coor_period=O2_coor;
                                        O_coor_period[1]-=y_length;
                                    }
                                    else if (O2_coor[1]-O1_coor[1]<-10){
                                        O_coor_period=O2_coor;
                                        O_coor_period[1]+=y_length;
                                    }

                                    if (judge_H_Bond(O1_coor, H1_coor, O_coor_period)) {
                                        connection_data[wat1_i][0]=O_WAT_id[wat2_i];

                                    }
                                    if (judge_H_Bond(O1_coor, H2_coor, O_coor_period)) {
                                        connection_data[wat1_i][1]=O_WAT_id[wat2_i];
                                    }
                                }
                                if (i==0 && j==0){

                                    if (O2_coor[1]-O1_coor[1]>10&&O2_coor[0]-O1_coor[0]>10){
                                        O_coor_period=O2_coor;
                                        O_coor_period[0]-=x_length;
                                        O_coor_period[1]-=y_length;
                                    }
                                    else  if (O2_coor[1]-O1_coor[1]<-10 && O2_coor[0]-O1_coor[0]<-10){
                                        O_coor_period=O2_coor;
                                        O_coor_period[0]+=x_length;
                                        O_coor_period[1]+=y_length;
                                    }

                                    if (judge_H_Bond(O1_coor, H1_coor, O_coor_period)) {
                                        connection_data[wat1_i][0]=O_WAT_id[wat2_i];
                                        //std::cout<<O_WAT_id[wat1_i]<<"   "<<O_WAT_id[wat2_i]<<"  ";
                                    }
                                    if (judge_H_Bond(O1_coor, H2_coor, O_coor_period)) {
                                        connection_data[wat1_i][1]=O_WAT_id[wat2_i];
                                        //std::cout<<O_WAT_id[wat1_i]<<"   "<<O_WAT_id[wat2_i]<<"  ";
                                    }
                                }
                            }
                        }
                    }
                }
            }
            for (index0 i = 0; i != O_WAT_id.size(); ++i){
                outfile<<O_WAT_id[i]<<std::setw(7)<<nc<<std::setw(7)<<frame<<std::setw(7)<<connection_data[i][0]<<std::setw(7)<<connection_data[i][1]<<std::endl;
            }
        }
        outfile.close();
    }
}

