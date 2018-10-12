//
// Created by utena on 17-11-21.
//
#include "Configurations&Functions.hpp"
#include "Rotation_Relaxation_HBNum_test.hpp"
#include "Diffusion_HBNum_test.hpp"
#include "Rotation_Relaxation&Diffusion_test.hpp"
#include "RotationalRelaxation&DiffusionByCluster.hpp"
#include "ExtractInterlayerTransportTimeSeriesData.hpp"
int main(int argc, char* argv[]){
    char dist[16];
    std::cout<<argv[1]<<std::endl;
    sprintf(dist,"%s",argv[1]);
    int temperature=std::atoi(argv[2]);
    int start_nc=std::atoi(argv[3]);
    int end_nc=std::atoi(argv[4]);
    //CalcHBondNetwork(dist, temperature, start_nc, end_nc);
    //CalcClusterByHBond(dist,temperature,start_nc,end_nc);
    //RCF_HBNum(dist,temperature,start_nc,end_nc);
    //Diffusion_HBNum(dist, temperature, start_nc, end_nc);
    //CalcRotationalRelaxationCorrelationFunction(dist, temperature, start_nc, end_nc);
    ExtractWaitingData(start_nc, end_nc);
    return 0;
}
