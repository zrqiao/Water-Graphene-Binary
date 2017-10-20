//
// Created by utena on 17-10-14.
//

#ifndef QZR_IDENTIFICATION_OF_GRID_WATER_H
#define QZR_IDENTIFICATION_OF_GRID_WATER_H


#include "../amber_netcdf.hpp"
#include "../amber_parm_1.0.hpp"
#include "../vector_calc.h"

typedef std::vector<double>::size_type index;

double calc_angle_z_axis(std::vector<double> ori){
    std::vector<double> z_axis={0,0,1};
    double angle;
    angle = vector_calc::vector_angle(z_axis,ori);
    return 180.0*acos(angle)/vector_calc::PI;
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