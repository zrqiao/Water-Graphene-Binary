//
// Created by utena on 17-10-14.
//

#ifndef QZR_IDENTIFICATION_OF_GRID_WATER_H
#define QZR_IDENTIFICATION_OF_GRID_WATER_H

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
#define x_length 99.4194366
#define y_length 96.5999970
#define dt 1

#define dr 0.01
#define C_NUM 3772
#define name_parm7 "../../../density_dis6a5_WAT1150.parm7"
//#define dist "7a5"
//#define temperature 320
//~ #define name_nc "water_ion_graphene_10a5"c
typedef std::vector<double>::size_type index0;



namespace dir_tool
{
    int mkpath(std::string s,mode_t mode=0755)
    {
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
int CreatDir(char *pDir)
{
    int i = 0;
    int iRet;
    int iLen;
    char* pszDir;

    if(NULL == pDir)
    {
        return 0;
    }

    pszDir = strdup(pDir);
    iLen = strlen(pszDir);

    // 创建中间目录
    for (i = 0;i < iLen;i ++)
    {
        if (pszDir[i] == '\\' || pszDir[i] == '/')
        {
            pszDir[i] = '\0';

            //如果不存在,创建
            iRet = ACCESS(pszDir,0);
            if (iRet != 0)
            {
                iRet = MKDIR(pszDir);
                if (iRet != 0)
                {
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

double calc_angle_z_axis(std::vector<double> ori) {
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



#endif //QZR_IDENTIFICATION_OF_GRID_WATER_H