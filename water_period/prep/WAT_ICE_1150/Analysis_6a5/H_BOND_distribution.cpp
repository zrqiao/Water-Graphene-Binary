#include "../amber_netcdf.hpp"
#include "../amber_parm_1.0.hpp"
#include "../vector_calc.h"
#include "Identification_of_grid_water.hpp"
//~ #include <string.h>
//~ #include <fstream>
//~ #include <iostream>
//~ #include <stdlib.h>
//~ #include <stdio.h>


#define start_nc 0
#define end_nc  20
#define select_boundary 3.55
#define hbond_cutoff_up 3.5
#define hbond_cutoff_down 2.0
#define hbond_cutoff_angle_up 180.0
#define hbond_cutoff_angle_down 135.0
#define period_boundary
#define x_length 99.4194366
#define y_length 96.5999970
#define dt 1
#define frame_segment 100
#define dr 0.01
#define name_parm7 "density_dis6a5_WAT1150.parm7"
#define C_NUM 3772
//~ #define name_nc "water_ion_graphene_10a5"

int main() {
	typedef std::vector<double>::size_type index;
	std::cout << "program to calculate the H-BOND and cluster distribution"<< "\n"<< std::endl;
    char HB_distribution_out_name[64];
    sprintf(HB_distribution_out_name, "Analysis_6a5/HB_distribution/310K/nc%d-%d_25ps.dat", start_nc,end_nc);
    std::ofstream outfile;
    outfile.open(HB_distribution_out_name);
    for (int nc = start_nc; nc <= end_nc; nc++) {
        char infile_name[64];
        sprintf(infile_name, "Analysis_6a5/HB_Network/310K/nc_%d.dat", nc);
        char name_nc[64];
        sprintf(name_nc, "310_0a25ps_TOT200ns/density_dis6a5_WAT1150_NVT_310_%d.nc", nc);
        amber_parm parm_name(name_parm7);
        std::vector <index> O_WAT_id = parm_name.id_by_type("OW");
        std::ifstream infile;
        
        int tot_frame=10000;
        for (int segment_id=0;segment_id<tot_frame/frame_segment;segment_id++) {
            int HB_Number_SUM[5]={0};

            //~ double total_hbond = 0;
            double total_wat =0;
            std::vector<std::vector<std::vector<double>>> nc_HB_connection_list(frame_segment,std::vector<std::vector<double>>(15000,std::vector<double>()));
            std::cout<<nc<<"  "<<segment_id<<std::endl;
            int start_frame= segment_id*frame_segment;
            int end_frame=(segment_id+1)*frame_segment;
            double num[5];
            infile.open(infile_name);
            while (!infile.fail()){
                infile>>num[0]>>num[1]>>num[2]>>num[3]>>num[4];
                if (num[2]>=start_frame&&num[2]<end_frame) {
                    double Atom0 = num[0];
                    double nc = num[1];
                    int frame = num[2]-start_frame;
                    //std::cout<<num[3]<<std::endl;
                    if (num[3] != -1) {
                        nc_HB_connection_list[frame][Atom0].push_back(num[3]);
                        nc_HB_connection_list[frame][num[3]].push_back(Atom0);
                    }
                    if (num[4] != -1) {
                        nc_HB_connection_list[frame][Atom0].push_back(num[4]);
                        nc_HB_connection_list[frame][num[4]].push_back(Atom0);
                    }
                }
            }
            infile.close();
            for (int framer = 0; framer < frame_segment; framer += dt) {
                for (int i = 0; i != O_WAT_id.size(); i++) {
                    int HB_num=nc_HB_connection_list[framer][O_WAT_id[i]].size();
                    HB_Number_SUM[HB_num]+=1;
                    total_wat++;
                }
            }
            /*for (int frame = 0; frame < frame_segment; frame += dt) {
                for (int i = 0; i != O_WAT_id.size(); i++) {
                    nc_HB_connection_list[frame][O_WAT_id[i]].clear();
                }
            }*/
            for (int HBNUM=0;HBNUM<5;HBNUM++){
                outfile<<HB_Number_SUM[HBNUM]/total_wat<<std::setw(15);
            }
            outfile<<std::endl;


        }

        //std::cout<<"finished";
        /*
        int Target_ID, frame, HB_ID1, HB_ID2;
        int x_points=ceil(x_length/select_boundary)+1;
        int y_points=ceil(y_length/select_boundary)+1;
        for (index frame = 0; frame < nc_data.frames_number(); frame += dt) {
            std::cout<<nc<<' '<<frame<<std::endl;
            std::vector < std::vector < std::vector < int > > > xypoints(x_points+1, std::vector < std::vector < int >> (y_points+1, std::vector<int>(0)));
            std::vector<std::vector<int>> check_connection(O_WAT_id.size(),std::vector<int>(O_WAT_id.size(),0));
            std::vector<std::vector<int>> connection_data(O_WAT_id.size(),std::vector<int>(2,-1));
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
                x    long long total_water_SUM;ypoints[0][j].insert(xypoints[0][j].end(), xypoints[x_points-1][j].begin(), xypoints[x_points-1][j].end());
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
            for (index i = 0; i != O_WAT_id.size(); ++i){
                outfile<<O_WAT_id[i]<<std::setw(7)<<nc<<std::setw(7)<<frame<<std::setw(7)<<connection_data[i][0]<<std::setw(7)<<connection_data[i][1]<<std::endl;
            }
        }
        outfile.close();
        */
    }
    outfile.close();
}


