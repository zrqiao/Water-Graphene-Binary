//
// Created by utena on 17-10-19.
//

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
#define frame_segment 1
#define dr 0.01
#define name_parm7 "density_dis6a5_WAT1150.parm7"
#define C_NUM 3772
//~ #define name_nc "water_ion_graphene_10a5"

int main() {
    typedef std::vector<double>::size_type index;
    std::cout << "program to calculate the H-BOND and cluster distribution"<< "\n"<< std::endl;
    char out_name[64];
    amber_parm parm_name(name_parm7);
    std::vector <index> O_WAT_id = parm_name.id_by_type("OW");
    long long cluster_size_distribution[O_WAT_id.size()]={0};
    for (int nc = start_nc; nc <= end_nc; nc++) {
        char infile_name[64];
        sprintf(infile_name, "Analysis_6a5/HB_Network/310K/nc_%d.dat", nc);
        char name_nc[64];
        sprintf(name_nc, "310_0a25ps_TOT200ns/density_dis6a5_WAT1150_NVT_310_%d.nc", nc);

        int WAT_ID_to_i[12000]={-1};
        for (int i=0; i<12000;i++){
            WAT_ID_to_i[i]=-1;
        }
        for (int i=0; i<O_WAT_id.size();i++){
            WAT_ID_to_i[O_WAT_id[i]]=i;
        }
        std::ifstream infile;
        int tot_frame=10000;

        sprintf(out_name, "Analysis_6a5/Cluster_Data/310K/nc_%d.dat", nc);
        std::ofstream outfile;
        outfile.open(out_name);
        double num[5];
        infile.open(infile_name);
        infile>>num[0]>>num[1]>>num[2]>>num[3]>>num[4];
        for (int frame_id=0;frame_id<tot_frame;frame_id++) {
            int HB_Number_SUM[5]={0};

            //~ double total_hbond = 0;
            double total_wat =0;
            std::vector<int> nc_HB_connection_list[O_WAT_id.size()]={std::vector<int>()};
            //std::cout<<nc<<"  "<<frame_id<<std::endl;
            int start_frame= frame_id;
            int end_frame=frame_id+1;

            while (!infile.fail()){
                int Atom0 = num[0];
                int nc = num[1];
                int frame = num[2]-start_frame;
                int O1=num[3];
                int O2=num[4];
                //std::cout<<num[3]<<std::endl;
                if (num[3] != -1) {
                    nc_HB_connection_list[WAT_ID_to_i[Atom0]].push_back(O1);
                    nc_HB_connection_list[WAT_ID_to_i[O1]].push_back(Atom0);
                }
                if (num[4] != -1) {
                    nc_HB_connection_list[WAT_ID_to_i[Atom0]].push_back(O2);
                    nc_HB_connection_list[WAT_ID_to_i[O2]].push_back(Atom0);
                }
                if (num[2]>=start_frame&&num[2]<end_frame) {
                    infile>>num[0]>>num[1]>>num[2]>>num[3]>>num[4];
                }
                if (num[2]>=end_frame){break;}
            }

            //统计团簇归类(搜索+并查集)
            std::vector<int> temp_path;
            std::vector<int> temp_succ_path;


            std::vector<int> Union_Find_children[O_WAT_id.size()]={std::vector<int>()};//存储的是O_WAT_id的index
            int Union_Find_father[O_WAT_id.size()]={0};
            for (int i = 0; i != O_WAT_id.size(); i++) {
                Union_Find_father[i]=i;
                Union_Find_children[i].push_back(i);
            }
            for (int i = 0; i != O_WAT_id.size(); i++) {
                //std::cout<<O_WAT_id[i]<<std::endl;
                temp_path.clear();
                temp_succ_path.clear();
                SearchGrid(O_WAT_id[i], O_WAT_id[i], temp_path, nc_HB_connection_list, temp_succ_path, WAT_ID_to_i);
                nc_HB_connection_list[i].clear();
                for (std::vector<int>::iterator iter=temp_succ_path.begin();iter!=temp_succ_path.end();iter++){
                    int conn_wat_i=WAT_ID_to_i[*iter];
                    if (Union_Find_father[i]!=Union_Find_father[conn_wat_i]){//如果不连通
                        int Insert_Target=Union_Find_father[i];//加入本组并查集
                        int Insert_Object=Union_Find_father[conn_wat_i];
                        for (std::vector<int>::iterator iter1=Union_Find_children[Insert_Object].begin();iter1!=Union_Find_children[Insert_Object].end();iter1++){
                            int object=*iter1;
                            Union_Find_father[object]=Insert_Target;//连接全部节点至根节点1
                            //std::cout<<O_WAT_id[object]<<std::setw(5);
                        }
                        //std::cout<<std::endl;
                        Union_Find_children[Insert_Target].insert(Union_Find_children[Insert_Target].end(),Union_Find_children[Insert_Object].begin(),Union_Find_children[Insert_Object].end());//整体重定向
                        Union_Find_children[Insert_Object].clear();
                        //Union_Find_children[Insert_Object].push_back(Insert_Object);
                    }
                }
            }
            outfile<<"#nc"<<nc<<" #frame"<<frame_id<<std::endl;
            std::cout<<"#nc"<<nc<<" #frame"<<frame_id<<std::endl;
            for (int i = 0; i != O_WAT_id.size(); i++) {
                //std::cout<<i<<std::setw(7);
                outfile<<i<<std::setw(7);
                for (std::vector<int>::iterator iter1=Union_Find_children[i].begin();iter1!=Union_Find_children[i].end();iter1++){
                    int object=*iter1;
                    //std::cout<<O_WAT_id[object]<<std::setw(7);
                    outfile<<O_WAT_id[object]<<std::setw(7);
                }
                cluster_size_distribution[Union_Find_children[i].size()]++;
                outfile<<std::endl;
                //std::cout<<std::endl;
            }
            std::ofstream distribution_out;
            distribution_out.open("Analysis_6a5/Cluster_Data/310K/cluster_size_distribution_nc0-20.dat");
            for (int i = 0; i != O_WAT_id.size(); i++) {
                distribution_out<<i<<std::setw(15)<<cluster_size_distribution[i]<<std::endl;
            }
            distribution_out.close();

        }
        infile.close();
        outfile.close();
    }
}


