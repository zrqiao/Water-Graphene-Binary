//
// Created by utena on 17-10-14.
//

#ifndef QZR_IDENTIFICATION_OF_GRID_WATER_H
#define QZR_IDENTIFICATION_OF_GRID_WATER_H


#include "amber_netcdf.hpp"
#include "amber_parm_1.0.hpp"
#include "vector_calc.h"

typedef std::vector<double>::size_type index;

double calc_angle_z_axis(std::vector<double> ori){
    std::vector<double> z_axis={0,0,1};
    double angle;
    angle = vector_calc::vector_angle(z_axis,ori);
    return 180.0*acos(angle)/vector_calc::PI;
}
std::vector<std::vector<std::vector<double> > > generate_HB_connection_list(char infile_name, int tot_frame){
    std::vector< std::vector< std::vector <double > > > new_connection_list(tot_frame,std::vector<std::vector<double>>(10000,new std::vector<double>));
    double num[5];
    infile=std::ifstream.open(infile_name);

    while (!infile.fail()){
        infile>>num[0]>>num[1]>>num[2]>>num[3]>>num[4];
        double Atom0=num[0];
        double nc=num[1];
        double frame=num[2];
        if (num[3]!=-1){
            new_connection_list[frame][Atom0].push_back(num[3]);
            new_connection_list[frame][num[3]].push_back(Atom0);
        }
        if (num[4]!=-1){
            new_connection_list[frame][Atom0].push_back(num[4]);
            new_connection_list[frame][num[4]].push_back(Atom0);
        }
    }
    return new_connection_list;
}

int SearchGrid(double Atom0,double Target,std::vector<double> Path,double frame,std::vector<std::vector<std::vector>> connection_list,std::vector<double> &success_path){
    for (int i=0;i<connection_list[frame][Atom0].size();i++){
        if (connection_list[frame][Atom0][i]==Target){
            if (Path.size()==4){
                for (j=0;j<3;j++){
                    success_path.push_back(Path[j]);
                }
            }
            return Path;
        } else if (Path.size()<4){
            std::vector<double> newPath=Path;
            newPath.push_back(connection_list[frame][Atom0][i]);
            SearchGrid(connection_list[frame][Atom0][i],Target,newPath,frame,connection_list);
        }
    }
    return success_path;
}



#endif //QZR_IDENTIFICATION_OF_GRID_WATER_H