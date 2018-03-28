//
// Created by utena on 17-10-19.
//

#include "Identification_of_grid_water.hpp"
//~ #include <string.h>
//~ #include <fstream>
//~ #include <iostream>
//~ #include <stdlib.h>
//~ #include <stdio.h>

int Calc_Cluster(char* dist,int temperature,int start_nc,int end_nc) {
    CreatDir("Analysis/Cluster_Data");
    typedef std::vector<double>::size_type index;
    std::cout << "program to calculate the H-BOND and cluster distribution"<< "\n"<< std::endl;
    char out_name[64];
    amber_parm parm_name(name_parm7);
    std::vector <index0> O_WAT_id = parm_name.id_by_type("OW");
    int num[5];
    int WAT_ID_to_i[12000];
    for (int i=0; i<12000;i++) WAT_ID_to_i[i]=-1;
    for (int i=0; i<O_WAT_id.size();i++) WAT_ID_to_i[O_WAT_id[i]]=i;
    for (int nc = start_nc; nc <= end_nc; nc++) {
        long double cluster_size_distribution[O_WAT_id.size()]={0};
        char infile_name[64];
        sprintf(infile_name, "Analysis/HB_Network/nc_%02d.dat", nc);
        std::ifstream infile;
        infile.open(infile_name);

        sprintf(out_name, "Analysis/Cluster_Data/nc_%02d.dat", nc);
        std::ofstream outfile;
        outfile.open(out_name);
        std::ofstream max_cluster_out;
        char max_out_name[64];
        sprintf(max_out_name,"Analysis/Cluster_Data/max_cluster_size_nc%02d.dat",nc);
        max_cluster_out.open(max_out_name);
        char name_nc[64];
        sprintf(name_nc, "density_dis%s_WAT1150_NVT_%d_%d.nc",dist,temperature,nc);
        nctraj nc_data(name_nc);
        int tot_frame=nc_data.frames_number();
        long long tot_wat=tot_frame*O_WAT_id.size();
        infile>>num[0]>>num[1]>>num[2]>>num[3]>>num[4];
        for (int frame_id=0;frame_id<tot_frame;frame_id++) {
            //int HB_Number_SUM[5]={0};

            //~ double total_hbond = 0;
            std::vector<int> nc_HB_connection_list[O_WAT_id.size()]={std::vector<int>()};
            //std::cout<<nc<<"  "<<frame_id<<std::endl;
            int Atom0,nc,frame,O1,O2;
            while (!infile.fail()){
                //std::cout<<num[0]<<' '<<num[1]<<' '<<num[2]<<" "<<num[3]<<" "<<num[4]<<std::endl;
                Atom0 = num[0];
                nc = num[1];
                frame = num[2];
                O1=num[3];
                O2=num[4];
                //std::cout<<Atom0<<' '<<nc<<' '<<num[2]<<" "<<num[3]<<" "<<num[4]<<std::endl;

                //std::cout<<num[3]<<std::endl;

                if (O2 != -1) {//condition- Take care!
                    nc_HB_connection_list[WAT_ID_to_i[Atom0]].push_back(O2);
                    //std::cout<<num[0]<<' '<<num[1]<<' '<<num[2]<<" "<<num[3]<<" "<<num[4]<<std::endl;

                    nc_HB_connection_list[WAT_ID_to_i[O2]].push_back(Atom0);
                    //std::cout<<num[0]<<' '<<num[1]<<' '<<num[2]<<" "<<num[3]<<" "<<num[4]<<std::endl;
                }
                if (O1 != -1) {
                    nc_HB_connection_list[WAT_ID_to_i[Atom0]].push_back(O1);
                    nc_HB_connection_list[WAT_ID_to_i[O1]].push_back(Atom0);
                }

                if (frame==frame_id) {

                    infile>>num[0]>>num[1]>>num[2]>>num[3]>>num[4];

                }
                if (frame>frame_id){break;}
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
            int cluster_size;
            int max_cluster_size=0;
            int max2_cluster_size=0;
            int max3_cluster_size=0;
            for (int i = 0; i != O_WAT_id.size(); i++) {
                //std::cout<<i<<std::setw(7);
                outfile<<O_WAT_id[i]<<std::setw(7);
                /*for (std::vector<int>::iterator iter1=Union_Find_children[i].begin();iter1!=Union_Find_children[i].end();iter1++){
                    int object=*iter1;
                    //std::cout<<O_WAT_id[object]<<std::setw(7);
                    outfile<<O_WAT_id[object]<<std::setw(7);
                }*/
                outfile<<Union_Find_father[i];
                cluster_size=Union_Find_children[i].size();
                cluster_size_distribution[cluster_size]++;
                if (cluster_size>max_cluster_size){
                    max3_cluster_size=max2_cluster_size;
                    max2_cluster_size=max_cluster_size;
                    max_cluster_size=cluster_size;
                }
                outfile<<std::endl;
                //std::cout<<std::endl;
            }
            max_cluster_out<<nc*tot_frame+frame_id<<std::setw(7)<<max_cluster_size<<std::setw(7)<<max2_cluster_size<<std::setw(7)<<max3_cluster_size<<std::endl;
        }
        infile.close();
        outfile.close();

        std::ofstream distribution_out;
        sprintf(out_name, "Analysis/Cluster_Data/cluster_size_distribution_nc%02d.dat", nc);
        distribution_out.open(out_name);
        for (int i = 1; i != O_WAT_id.size(); i++) {
            distribution_out<<i<<std::setw(15)<<cluster_size_distribution[i]/tot_wat<<std::endl;
        }
        distribution_out.close();
        max_cluster_out.close();
    }
}


