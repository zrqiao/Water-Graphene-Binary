//
// Created by utena on 17-6-14.
//
#include "enthaply/amber_netcdf.hpp"
#include "enthaply/amber_parm_1.0.hpp"
#include<cmath>
#include<cstdlib>
#include<cstring>
#include<cfloat>
#include <algorithm>
#define dx 0.2
int main(){
    std::ifstream infile;
    std::ofstream outfile;
    infile.open("cavity_size_change");
    outfile.open("size_change_distribution");
    std::vector<int> size_change_distribution(4000,0);
    while (!infile.fail()) {
        double num[4];
        double cavity_size_change;
        infile >> num[0] >> num[1] >> num[2] >> num[3];
        cavity_size_change=num[3];
        std::cout<<floor(cavity_size_change / pow(dx,2)+3000)<<std::endl;
        size_change_distribution[floor(cavity_size_change / pow(dx,2)+3000)]++;
    }
    for (int i=0;i<4000;i++) {
        std::cout<<(i-3000)*0.04<<std::setw(7)<<size_change_distribution[i]<<std::endl;
    }
    infile.close();
    outfile.close();
}
