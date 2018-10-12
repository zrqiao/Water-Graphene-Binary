//
// Created by utena on 17-10-14.
//

#ifndef QZR_CONFIGURATIONS
#define QZR_CONFIGURATIONS

#include "amber_netcdf.hpp"
#include "amber_parm_1.0.hpp"
#include "vector_calc.h"
#include "unistd.h"
#include <stdarg.h>
    #include <sys/stat.h>
#define ACCESS access
    #define MKDIR(a) mkdir((a),0755)
#include <dirent.h>
#include <stdlib.h>
#include <string.h>
//#define start_nc 0
//#define end_nc  11
#define select_boundary 3.55
#define hbond_cutoff_up 3.5
#define hbond_cutoff_down 2.0
#define hbond_cutoff_angle_up 180.0
#define hbond_cutoff_angle_down 135.0
#define period_boundary
#define x_length  87.0363607
#define y_length  85.9673587
#define dt 4 //16fs
#define deviation_y_center 30
#define deviation_x_center 30
#define dr 0.01
#define C_NUM 1400
#define parm7_prefix "density_dis9a5"
#define name_parm7 "nc/density_dis9a5.parm7"
#define nc_path "nc"
#define width_transition_path 1.2
//#define dist "7a5"
//#define temperature 320
//~ #define name_nc "water_ion_graphene_10a5"c
typedef std::vector<double>::size_type index0;
int ReduceFrame(std::vector<std::vector<int>> frame_sta, std::vector<std::vector<int>> &condensed_frames);
int PickFrameByCluster(int frame_zeropoint, std::vector<int> Cluster_Size_by_WATID_frame,
                       std::vector<std::vector<int>> &condensed_state_frames);
int PickFrameByLayer(int start_nc, int end_nc, const std::vector<double>& boundary_vector,
                     index0 o_wat_id, std::vector<std::vector<int>> &condensed_state_frames);

namespace dir_tool {
    int mkpath(std::string s,mode_t mode=0755) {
        size_t pre=0,pos;
        std::string dir;
        int mdret;

        if(s[s.size()-1]!='/'){
            // force trailing / so we can handle everything in loop
            s+='/';
        }

        while((pos=s.find_first_of('/',pre))!=std::string::npos){
            dir=s.substr(0,pos++);
            pre=pos;
            if(dir.size()==0) continue; // if leading / first time is 0 length
            if((mdret=::mkdir(dir.c_str(),mode)) && errno!=EEXIST){
                return mdret;
            }
        }
        return mdret;
    }
/*
 * 注：使用时，light::mkdir, 易与::mkdir 混淆，固去掉了
int mkdir(std::string s,mode_t mode)
{
    light::mkpath(s, mode);
}
*/
}
int CreatDir(char *pDir) {
    int i = 0;
    int iRet;
    int iLen;
    char* pszDir;
    if(NULL == pDir) {
        return 0;
    }
    pszDir = strdup(pDir);
    iLen = strlen(pszDir);
    // 创建中间目录
    for (i = 0;i < iLen;i ++) {
        if (pszDir[i] == '\\' || pszDir[i] == '/')
        {
            pszDir[i] = '\0';
            //如果不存在,创建
            iRet = ACCESS(pszDir,0);
            if (iRet != 0) {
                iRet = MKDIR(pszDir);
                if (iRet != 0) {
                    return -1;
                }
            }
            //支持linux,将所有\换成/
            pszDir[i] = '/';
        }
    }
    iRet = MKDIR(pszDir);
    free(pszDir);
    return iRet;
}

double CalcAngleZAxis(std::vector<double> ori) {
    std::vector<double> z_axis = {0, 0, 1};
    double angle;
    angle = vector_calc::vector_angle(z_axis, ori);
    return 180.0 * acos(angle) / vector_calc::PI;
}

double generate_HB_connection_list(std::vector<std::vector<std::vector<double>>> &new_connection_list, char* infile_name, int start_frame,int end_frame){
    double num[5];
    std::ifstream infile;
    infile.open(infile_name);
    while (!infile.fail()){
        infile>>num[0]>>num[1]>>num[2]>>num[3]>>num[4];
        if (num[2]>=start_frame&&num[2]<end_frame) {
            double Atom0 = num[0];
            double nc = num[1];
            double frame = num[2]-start_frame;
            if (num[3] != -1) {
                new_connection_list[frame][Atom0].push_back(num[3]);
                new_connection_list[frame][num[3]].push_back(Atom0);
            }
            if (num[4] != -1) {
                new_connection_list[frame][Atom0].push_back(num[4]);
                new_connection_list[frame][num[4]].push_back(Atom0);
            }
        }
    }
    infile.close();
    return 1;
}

bool In(int Target,std::vector<int> array){
    for (int i=0;i<array.size();i++){
        if (Target==array[i]){
            return true;
        }
    }
    return false;
}

std::vector<int> SearchGrid(int Atom0, int Target, std::vector<int> Path, std::vector<int>* connection_list,
                            std::vector<int> &success_path, int *WAT_ID_to_i) {
    for (int i=0;i<connection_list[WAT_ID_to_i[Atom0]].size();i++){
        if (connection_list[WAT_ID_to_i[Atom0]][i]==Target){
            if (Path.size()==3){
                for (int j=0;j<3;j++){
                    success_path.push_back(Path[j]);
                }
            }
        } else if (Path.size()<3){
            std::vector<int> newPath=Path;
            newPath.push_back(connection_list[WAT_ID_to_i[Atom0]][i]);
            int newtarget=connection_list[WAT_ID_to_i[Atom0]][i];
            if (!In(newtarget,Path)) {
                SearchGrid(newtarget, Target, newPath, connection_list, success_path,WAT_ID_to_i);
            }
        }
    }
    return success_path;
}

inline int JudgeLayer(const std::vector<double>& o_coor, double bilayer_center){

    if(o_coor[2] > bilayer_center + width_transition_path){
        return 1;
    }
    else if(o_coor[2] < bilayer_center - width_transition_path){
        return -1;
    }
    else{
        return 0;
    }
}

inline int JudgeHBond(std::vector<double> O1_coor, std::vector<double> H_coor, std::vector<double> O2_coor){
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

std::vector<double> InitializeBoundary(int start_nc, int end_nc){
    char name_nc[64];
    double C_x_coor_sum = 0, C_y_coor_sum = 0;
    double C_x_coor_center, C_y_coor_center;
    double C_z_coor_sum_1 = 0, C_z_coor_sum_2 = 0;
    double C_z_coor_average_1, C_z_coor_average_2;
    double Z_UP,Z_DOWN,Y_UP,Y_DOWN,X_UP,X_DOWN;
    for (int nc=start_nc;nc<=end_nc;nc++) {
        sprintf(name_nc, "%s/%s_%d.nc",nc_path, parm7_prefix, nc);
        nctraj data_nc(name_nc);
        for (index0 C_index = 0; C_index != C_NUM; ++C_index) {
            C_z_coor_sum_1 += data_nc.atom_coordinate(0, C_index)[2];
            C_z_coor_sum_2 += data_nc.atom_coordinate(0, (C_NUM + C_index))[2];
            C_y_coor_sum += data_nc.atom_coordinate(0, C_index)[1];
            C_x_coor_sum += data_nc.atom_coordinate(0, C_index)[0];
        }
    }
    C_z_coor_average_1 = C_z_coor_sum_1/(C_NUM*(end_nc-start_nc+1));
    C_z_coor_average_2 = C_z_coor_sum_2/(C_NUM*(end_nc-start_nc+1));
    C_y_coor_center = C_y_coor_sum/(C_NUM*(end_nc-start_nc+1));
    C_x_coor_center = C_x_coor_sum/(C_NUM*(end_nc-start_nc+1));

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

    std::vector<double> boundary_vector={X_UP,X_DOWN,Y_UP,Y_DOWN,Z_UP,Z_DOWN};
    return boundary_vector;
}

bool IsInGraphene(std::vector<double> boundary_vector, std::vector<double> coor){
    double Z_UP,Z_DOWN,Y_UP,Y_DOWN,X_UP,X_DOWN;
    X_UP=boundary_vector[0];
    X_DOWN=boundary_vector[1];
    Y_UP=boundary_vector[2];
    Y_DOWN=boundary_vector[3];
    Z_UP=boundary_vector[4];
    Z_DOWN=boundary_vector[5];
    if(coor[0]>X_DOWN&&coor[0]<X_UP&&coor[1]>Y_DOWN&&coor[1]<Y_UP&&coor[2]>Z_DOWN&&coor[2]<Z_UP){
        return true;
    } else{
        return false;
    }
}

int CalcHBondNetwork(char *dist, int temperature, int start_nc, int end_nc) {
    CreatDir("Analysis/HB_Network");
    std::vector<double> OH1_vec(3, 0), OH2_vec(3, 0), Dip_vec(3, 0), Nor_vec(3, 0), O1H1_vec(3, 0), O2H2_vec(3, 0);
    std::vector <index0> O_WAT_IN_C_id, select_wat_id;
    std::vector<double> O_coor, O1_coor, O2_coor, H1O_coor, H2O_coor, H1_coor, H2_coor;
    double angleH1, angleH2, angleO1H1, angleO2H2, angleDip, angleNor;
    double Z_UP,Z_DOWN,Y_UP,Y_DOWN,X_UP,X_DOWN;
    amber_parm parm_name(name_parm7);
    std::vector <index0> O_WAT_id = parm_name.id_by_type("OW");
    std::vector<double> boundary_vector=InitializeBoundary(start_nc, end_nc);
    for (int nc = start_nc; nc <= end_nc; nc++) {
        char out_name[64];
        sprintf(out_name, "Analysis/HB_Network/nc_%02d.dat",nc);
        std::ofstream outfile;
        outfile.open(out_name);
        char name_nc[64];
        sprintf(name_nc, "%s/density_dis%s_%d.nc", nc_path, dist, nc);
        nctraj nc_data(name_nc);
        int tot_frame=nc_data.frames_number();
        int watnum=O_WAT_id.size();
        std::cout<<O_WAT_id.size();
        int *check_connection[watnum];
        for (int i=0;i<watnum;i++) {
            check_connection[i]=new int[watnum];
        }
        int is_in_C[watnum];
        int Target_ID, frame, HB_ID1, HB_ID2;
        int x_points=ceil(x_length/select_boundary)+1;
        int y_points=ceil(y_length/select_boundary)+1;
        for (index0 frame = 0; frame < nc_data.frames_number(); frame += dt) {

            std::vector < int > *xypoints[x_points+1];
            for (int i=0;i<x_points+1;i++){
                xypoints[i]=new std::vector<int>[y_points+1];
            }
            std::cout<<nc<<' '<<frame<<std::endl;
            /*
            for (int i=0;i<x_points+1;i++){
                for (int j=0;j<y_points+1;j++){
                    xypoints[i][j]=new std::vector<int>;
            }
            */

            int connection_data[watnum][2];
            for (int i=0;i<watnum;i++){
                connection_data[i][0]=-1;
                connection_data[i][1]=-1;
                for (int j=0;j<watnum;j++){
                    check_connection[i][j]=0;
                }
            }
            //std::cout<<'checkpoint';
            if (frame % 100 == 0) std::cout << frame << std::endl;

            for (int i = 0; i != O_WAT_id.size(); i++) {//打格点
                //std::cout<<i<<std::endl;
                O_coor = nc_data.atom_coordinate(frame, O_WAT_id[i]);
                if(IsInGraphene(boundary_vector, O_coor)){
                    is_in_C[i]=1;
                    int x_location = floor((O_coor[0]-X_DOWN) / select_boundary)+1;
                    int y_location = floor((O_coor[1]-Y_DOWN) / select_boundary)+1;
                    xypoints[x_location][y_location].push_back(i);
                    xypoints[x_location + 1][y_location].push_back(i);
                    xypoints[x_location][y_location + 1].push_back(i);
                    xypoints[x_location + 1][y_location + 1].push_back(i);
                } else{
                    is_in_C[i]=0;
                }
            }
            /*
            for (int i = 0; i < x_points; i++) {
                xypoints[i][0].insert(xypoints[i][0].end(), xypoints[i][y_points-1].begin(), xypoints[i][y_points-1].end());
            }
            for (int j = 0; j < y_points; j++) {
                xypoints[0][j].insert(xypoints[0][j].end(), xypoints[x_points-1][j].begin(), xypoints[x_points-1][j].end());
            }
            xypoints[0][0].insert(xypoints[0][0].end(), xypoints[x_points-1][y_points-1].begin(), xypoints[x_points-1][y_points-1].end());
             */
            for (int i = 0; i < x_points; i++) {
                for (int j = 0; j < y_points; j++) {
                    for (std::vector<int>::iterator iter1 = xypoints[i][j].begin(); iter1 != xypoints[i][j].end(); iter1++) {
                        int wat1_i = *iter1;
                        for (std::vector<int>::iterator iter2 = xypoints[i][j].begin();
                             iter2 != xypoints[i][j].end(); iter2++) {
                            int wat2_i = *iter2;
                            if (check_connection[wat1_i][wat2_i]==0) {//避免重复计算
                                check_connection[wat1_i][wat2_i]=1;
                                O1_coor = nc_data.atom_coordinate(frame, O_WAT_id[wat1_i]);
                                O2_coor = nc_data.atom_coordinate(frame, O_WAT_id[wat2_i]);
                                H1_coor = nc_data.atom_coordinate(frame, O_WAT_id[wat1_i] + 1);
                                H2_coor = nc_data.atom_coordinate(frame, O_WAT_id[wat1_i] + 2);
                                if (JudgeHBond(O1_coor, H1_coor, O2_coor)) {
                                    connection_data[wat1_i][0]=O_WAT_id[wat2_i];
                                }
                                if (JudgeHBond(O1_coor, H2_coor, O2_coor)) {
                                    connection_data[wat1_i][1]=O_WAT_id[wat2_i];
                                }
                                std::vector<double> O_coor_period(3,0);
                                /*
                                if (i==0){//NOTE:This subroutine is for periodic boundary systems.
                                    /*
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
                                    /*
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
                                    /*
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
                                    if (JudgeHBond(O1_coor, H2_coor, O_coor_period)) {
                                        connection_data[wat1_i][1]=O_WAT_id[wat2_i];
                                        //std::cout<<O_WAT_id[wat1_i]<<"   "<<O_WAT_id[wat2_i]<<"  ";
                                    }
                                }
                                */
                            }
                        }
                    }
                }
            }
            for (index0 i = 0; i != O_WAT_id.size(); ++i){
                if (is_in_C[i]) outfile<<O_WAT_id[i]<<std::setw(7)<<nc<<std::setw(7)<<frame<<std::setw(7)<<connection_data[i][0]<<std::setw(7)<<connection_data[i][1]<<std::endl;
            }
        }
        outfile.close();
    }
}

int CalcClusterByHBond(char *dist, int temperature, int start_nc, int end_nc) {
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

int PickFrameByCluster(int frame_zeropoint, std::vector<int> Cluster_Size_by_WATID_frame,
                       std::vector<std::vector<int>> &condensed_state_frames) {
    std::vector<std::vector<int>> frame_in_graphene;
    for(int frame =0; frame < Cluster_Size_by_WATID_frame.size()-dt; frame+=dt)
    {
        //o_coor = nc_data.atom_coordinate(frame, o_wat_id);
        if(true) {//核心状态判断
            std::vector<int> temp;
            temp.push_back(frame + frame_zeropoint);//校正時間起點
            if ((Cluster_Size_by_WATID_frame[frame]+Cluster_Size_by_WATID_frame[frame+3]+ Cluster_Size_by_WATID_frame[frame+4]+
                 Cluster_Size_by_WATID_frame[frame+1] + Cluster_Size_by_WATID_frame[frame+2] )/5 > 12) {
                //temp.push_back(floor((Cluster_Size_by_WATID_frame[o_wat_id_ref][frame]) / coarse_grain_cutoff));
                temp.push_back(2);
            }
            else if ((Cluster_Size_by_WATID_frame[frame]+Cluster_Size_by_WATID_frame[frame+3]+ Cluster_Size_by_WATID_frame[frame+4]+
                      Cluster_Size_by_WATID_frame[frame+1] + Cluster_Size_by_WATID_frame[frame+2] )/5 > 3) {
                //temp.push_back(floor((Cluster_Size_by_WATID_frame[o_wat_id_ref][frame]) / coarse_grain_cutoff));
                temp.push_back(1);
            }
            else{temp.push_back(0);}
            frame_in_graphene.push_back(temp);
        }
    }
    ReduceFrame(frame_in_graphene, condensed_state_frames);
    frame_in_graphene.clear();
    return 0;
}

int PickFrameByHBondConf(int start_nc, int end_nc, const std::vector<double>& boundary_vector,
        std::vector<std::vector<int>> HB_connection_by_WATID_frame,
        index0 o_wat_id, std::vector<std::vector<int>> &condensed_state_frames) {
    std::vector<std::vector<int>> frame_in_graphene;
    amber_parm parm_name(name_parm7);
    double bilayer_center=(boundary_vector[5]+boundary_vector[4])/2;
    index0 frame_in_total_nc = 0;

    for(int nc = start_nc; nc  != end_nc+1; ++nc) {
        char name_nc[64];
        sprintf(name_nc, "%s/%s_%d.nc",nc_path, parm7_prefix, nc);
        nctraj nc_data(name_nc);
        std::vector<double>  o_coor;
        int frame_this_nc = nc_data.frames_number();
        for(int frame =0; frame < frame_this_nc; frame+=dt) {
            o_coor = nc_data.atom_coordinate(frame, o_wat_id);
            if( IsInGraphene(boundary_vector, o_coor)) {//核心状态判断
                std::vector<int>  buff;
                buff.push_back(frame_in_total_nc+frame);
                int HB_id1=HB_connection_by_WATID_frame[frame][0];
                int HB_id2=HB_connection_by_WATID_frame[frame][1];
                if(bilayer_center - width_transition_path < o_coor[2]
                && o_coor[2] > bilayer_center + width_transition_path
                &&HB_id1>0&&HB_id2>0){
                    if(JudgeLayer(nc_data.atom_coordinate(frame, HB_id1), bilayer_center) *
                          JudgeLayer(nc_data.atom_coordinate(frame, HB_id2), bilayer_center)==-1){
                        buff.push_back(1);
                    }
                    else{
                        buff.push_back(0);
                    }
                }
                else{
                    buff.push_back(0);
                }
                frame_in_graphene.push_back(buff);
            }
        }
        frame_in_total_nc = frame_in_total_nc + frame_this_nc;
    }
    ReduceFrame(frame_in_graphene, condensed_state_frames);
    frame_in_graphene.clear();
    return 0;
}

int PickFrameByLayer(int start_nc, int end_nc, const std::vector<double>& boundary_vector,
                                                index0 o_wat_id, std::vector<std::vector<int>> &condensed_state_frames) {
    std::vector<std::vector<int>> frame_in_graphene;
    amber_parm parm_name(name_parm7);
    double bilayer_center=(boundary_vector[5]+boundary_vector[4])/2;
    index0 frame_in_total_nc = 0;
    for(int nc = start_nc; nc  != end_nc+1; ++nc) {
        char name_nc[64];
        sprintf(name_nc, "%s/%s_%d.nc",nc_path, parm7_prefix, nc);
        nctraj nc_data(name_nc);
        std::vector<double>  o_coor;
        int frame_this_nc = nc_data.frames_number();
        for(int frame =0; frame < frame_this_nc; frame+=dt) {
            o_coor = nc_data.atom_coordinate(frame, o_wat_id);
            if( IsInGraphene(boundary_vector, o_coor)) {//核心状态判断
                std::vector<int>  buff;
                buff.push_back(frame_in_total_nc+frame);
                buff.push_back(JudgeLayer(o_coor, bilayer_center));
                frame_in_graphene.push_back(buff);
            }
        }
        frame_in_total_nc = frame_in_total_nc + frame_this_nc;
    }
    ReduceFrame(frame_in_graphene, condensed_state_frames);
    frame_in_graphene.clear();
    return 0;
}

int ReduceFrame(std::vector<std::vector<int>> frame_sta, std::vector<std::vector<int>> &condensed_frames) {
    //typedef std::vector<double>::size_type index;
    std::vector<int> temp_condensed_frames;//0-label 1-startframe 2-framelength 3-endframe
    for(int i=0; i != (frame_sta.size()); ++i) {
        if (i==0 )	{
            temp_condensed_frames.push_back(frame_sta[i][1]);
            temp_condensed_frames.push_back(frame_sta[i][0]);
            temp_condensed_frames.push_back(dt);
        }
        else if(frame_sta[i][0]!=frame_sta[i-1][0]+dt||frame_sta[i][1]!=frame_sta[i-1][1]){
            temp_condensed_frames.push_back(frame_sta[i-1][0]);
            condensed_frames.push_back(temp_condensed_frames);

            temp_condensed_frames.clear();
            temp_condensed_frames.push_back(frame_sta[i][1]);
            temp_condensed_frames.push_back(frame_sta[i][0]);
            temp_condensed_frames.push_back(dt);
        }
        else if(frame_sta[i][1]==frame_sta[i-1][1]){
            temp_condensed_frames[2]+=dt;
        }
        else if (i==frame_sta.size()-1){
            temp_condensed_frames.push_back(frame_sta[i][0]);
            condensed_frames.push_back(temp_condensed_frames);
            temp_condensed_frames.clear();
        }
    }
    return 0;
}

#endif //QZR_CONFIGURATIONS