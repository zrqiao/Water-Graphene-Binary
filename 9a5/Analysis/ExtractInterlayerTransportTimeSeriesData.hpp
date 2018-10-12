
#include "amber_netcdf.hpp"
#include "amber_parm_1.0.hpp"
#include "vector_calc.h"
#include "mpi.h"
#include "Configurations&Functions.hpp"
//~ #define z_axis_modify  2
//#define start_nc 1
//#define end_nc  10
//~ #define relaxation_time 25
//#define name_parm7 "density_dis9a5.parm7"
#define jump_time 5000
#define max_transition_frames 10000
#define ntwx 4
#define ddt 25

//int WAT_NUM[4]={1806,1700,2443,2550};
//~ #define name_nc "water_ion_graphene_10a5"
std::vector<std::vector<int>>  calc_jump(std::vector<std::vector<int>> condensed_state_frames);
std::vector<std::vector<int>>  calc_waitingtime(std::vector<std::vector<int>> condensed_state_frames);
int calc_distribution(std::vector<std::vector<int>> timelag_stat, std::vector<double> &time_distribution_array, int window);

int calc_distribution(std::vector<std::vector<int>> timelag_stat,std::vector<double> &time_distribution_array);
int ExtractTPTData(int start_nc, int end_nc){
    CreatDir("Analysis/Transport");
    typedef std::vector<double>::size_type index;
    std::vector<int *> jump_coor;
    std::ofstream TPT_distribution_outfile;
    std::ofstream TPT_outfile;
    char name_output[64];
    sprintf(name_output,"Analysis/Transport/transition_path_WATid_start_finish_direction_nc_%d-%d",start_nc,end_nc);
    TPT_outfile.open(name_output);
    std::vector<int> jump_time_distribution(jump_time, 0);
    std::cout << "program to calculate the transition path time series data" << "\n" << std::endl;
    std::vector<double> transitionpath_time_distribution(max_transition_frames, 0);
    std::vector<std::vector<int>> jump_stat;
    std::vector<double> boundary_vector=InitializeBoundary(start_nc, end_nc);
    amber_parm parm_nam(name_parm7);
    std::vector<std::vector<int>> pickframe;
    std::vector<index> O_WAT_id = parm_nam.id_by_type("OW");
    std::cout << "total water:  " << O_WAT_id.size() << "\n" << std::endl;
    std::vector<int> wat_jump_count;
    std::vector<std::vector<int>> wat_jump_distribution;

    for (index i = 0; i != O_WAT_id.size(); ++i) {
        std::cout << "water calculate and id: " << std::setw(8) << i + 1 << std::setw(8) << O_WAT_id[i]
                  << std::endl;
        PickFrameByLayer(start_nc, end_nc, boundary_vector, O_WAT_id[i], pickframe);
        //~ std::cout<<"frame size: " << pickframe.size()<<std::endl;
        for (index j = 0; j != pickframe.size(); ++j)
            std::cout << pickframe[j][0] << "   " << pickframe[j][1] << "   " << pickframe[j][2] << "   "
                      << pickframe[j][3] << std::endl;
        if (pickframe.size() != 0) {
            jump_stat = calc_jump(pickframe);
            for (int j = 0; j < jump_stat.size(); j++) {
                /*
                if (jump_stat[j][3] == -1) {
                }
                if (jump_stat[j][3] == 1) {
                }
                 */
                TPT_outfile << O_WAT_id[i] << std::setw(10) << jump_stat[j][0] << std::setw(10) << jump_stat[j][1]
                         << std::setw(10) << jump_stat[j][3] << std::endl;
            }
            calc_distribution(jump_stat, transitionpath_time_distribution, dt);
            jump_stat.clear();
        }
        pickframe.clear();
    }

    int frame_in_total_nc = 0;
    for (int nc = start_nc; nc != end_nc + 1; ++nc) {
        char name_nc[64];
        sprintf(name_nc, "%s/%s_%d.nc",nc_path, parm7_prefix, nc);
        amber_parm parm_name(name_parm7);
        nctraj nc_data(name_nc);
        std::vector<double> o_coor;
        int frame_this_nc = nc_data.frames_number();
        frame_in_total_nc = frame_in_total_nc + frame_this_nc;
    }
    sprintf(name_output,"Analysis/transition_path_time_distribution_nc_%d-%d",start_nc,end_nc);
    TPT_distribution_outfile.open(name_output);
    for (index i = 0; i < transitionpath_time_distribution.size(); i += 1) {
        //~ TPT_distribution_outfile <<std::setw(15) <<i<<std::setw(15)<<count_all_water[i] << std::endl;
        TPT_distribution_outfile << i * ntwx * dt << std::setw(15)
                << transitionpath_time_distribution[i] * 1000 / (frame_in_total_nc * ntwx) << std::setw(15)
                << std::endl;
    }
    TPT_distribution_outfile.close();
    TPT_outfile.close();
    return 0;
}

int ExtractWaitingData(int start_nc, int end_nc){
    typedef std::vector<double>::size_type index;
    std::vector<int *> jump_coor;
    std::ofstream waiting_time_distribution_outfile;
    std::ofstream Waiting_outfile;
    char name_output[64];
    CreatDir("Analysis/Transport");
    sprintf(name_output,"Analysis/Transport/waiting_WATid_start_finish_direction_nc_%d-%d",start_nc,end_nc);
    Waiting_outfile.open(name_output);
    std::vector<int> jump_time_distribution(jump_time, 0);
    std::cout << "program to calculate the transition path time series data" << "\n" << std::endl;
    std::vector<double> waiting_time_distribution(max_transition_frames, 0);
    std::vector<std::vector<int>> wait_stat;
    std::vector<double> boundary_vector=InitializeBoundary(start_nc, end_nc);
    amber_parm parm_nam(name_parm7);
    std::vector<std::vector<int>> pickframe;
    std::vector<index> O_WAT_id = parm_nam.id_by_type("OW");
    std::cout << "total water:  " << O_WAT_id.size() << "\n" << std::endl;
    std::vector<int> wat_jump_count;
    std::vector<std::vector<int>> wat_jump_distribution;
    for (index i = 0; i != O_WAT_id.size(); ++i) {
        std::cout << "water calculate and id: " << std::setw(8) << i + 1 << std::setw(8) << O_WAT_id[i]
                  << std::endl;
        PickFrameByLayer(start_nc, end_nc, boundary_vector, O_WAT_id[i], pickframe);
        std::cout<<"frame size: " << pickframe.size()<<std::endl;
        for (index j = 0; j != pickframe.size(); ++j)
            std::cout << pickframe[j][0] << "   " << pickframe[j][1] << "   " << pickframe[j][2] << "   "
                      << pickframe[j][3] << std::endl;
        if (pickframe.size() != 0) {
            wait_stat = calc_waitingtime(pickframe);
            for (int j = 0; j < wait_stat.size(); j++) {
                Waiting_outfile << O_WAT_id[i] << std::setw(10) << wait_stat[j][0] << std::setw(10) << wait_stat[j][1]
                            << std::setw(10) << wait_stat[j][3] << std::endl;
            }
            calc_distribution(wait_stat, waiting_time_distribution, ddt);
            wait_stat.clear();
        }
        pickframe.clear();
    }
    int frame_in_total_nc = 0;
    for (int nc = start_nc; nc != end_nc + 1; ++nc) {
        char name_nc[64];
        sprintf(name_nc, "%s/%s_%d.nc",nc_path, parm7_prefix, nc);
        amber_parm parm_name(name_parm7);
        nctraj nc_data(name_nc);
        std::vector<double> o_coor;
        int frame_this_nc = nc_data.frames_number();
        frame_in_total_nc = frame_in_total_nc + frame_this_nc;
    }
    sprintf(name_output,"Analysis/waiting_time_distribution_nc_%d-%d",start_nc,end_nc);
    waiting_time_distribution_outfile.open(name_output);
    for (index i = 0; i < waiting_time_distribution.size(); i += 1) {
        //~ waiting_time_distribution_outfile <<std::setw(15) <<i<<std::setw(15)<<count_all_water[i] << std::endl;
        waiting_time_distribution_outfile << i * ntwx * ddt << std::setw(15)
                                 << waiting_time_distribution[i] * (frame_in_total_nc * ntwx) << std::setw(15)
                                 << std::endl;
    }
    waiting_time_distribution_outfile.close();
    Waiting_outfile.close();
    return 0;
}

/*
std::vector<std::vector<int>>  Reduce(std::vector<std::vector<int>> frame_sta)
{
    typedef std::vector<double>::size_type index;
    std::vector<std::vector<int>> condensed_jump_frames;
    std::vector<int> temp_condensed_frames;//0-layer 1-startframe 2-framelength 3-endframe
    for(index i=0; i != (frame_sta.size()); ++i)
    {
        if (i==0 )	{
            temp_condensed_frames.push_back(frame_sta[i][1]);
            temp_condensed_frames.push_back(frame_sta[i][0]);
            temp_condensed_frames.push_back(dt);
        }
        else if(frame_sta[i][0]!=frame_sta[i-1][0]+dt||frame_sta[i][1]!=frame_sta[i-1][1]){
            temp_condensed_frames.push_back(frame_sta[i-1][0]);
            condensed_jump_frames.push_back(temp_condensed_frames);

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
            condensed_jump_frames.push_back(temp_condensed_frames);
            temp_condensed_frames.clear();
        }
    }

    return condensed_jump_frames;
}
 */


std::vector<std::vector<int>>  calc_jump( std::vector<std::vector<int>> condensed_state_frames){
    std::vector<std::vector<int>> jump_stat;
    std::vector<int> temp_jump;// 0-startframe 1-endframe 2-transitiontime 3-direction
    for (int i=1;i<condensed_state_frames.size()-1;i++){
        if (condensed_state_frames[i][0]==0){
            if(condensed_state_frames[i-1][3]+dt==condensed_state_frames[i][1]&&condensed_state_frames[i][3]+dt==condensed_state_frames[i+1][1]){//必须完全连接
                if (condensed_state_frames[i-1][0]*condensed_state_frames[i+1][0]==-1){
                    temp_jump.push_back(condensed_state_frames[i][1]);
                    temp_jump.push_back(condensed_state_frames[i][3]);
                    temp_jump.push_back(condensed_state_frames[i][2]);
                    temp_jump.push_back(condensed_state_frames[i+1][0]);
                    jump_stat.push_back(temp_jump);
                    temp_jump.clear();
                }
            }
        }
    }
    return jump_stat;
}

std::vector<std::vector<int>>  calc_waitingtime( std::vector<std::vector<int>> condensed_state_frames){
    std::vector<std::vector<int>> waiting_stat;
    std::vector<int> temp_waiting;// 0-startframe 1-endframe 2-transitiontime 3-direction
    int last_success_frame= -1;
    for (int i=1;i<condensed_state_frames.size()-1;i++){
        if (condensed_state_frames[i-1][3]+dt==condensed_state_frames[i][1]&&condensed_state_frames[i][3]+dt==condensed_state_frames[i+1][1]){
            if(condensed_state_frames[i][0]==0){//必须完全连接
                if (condensed_state_frames[i-1][0]*condensed_state_frames[i+1][0]==-1){
                    if (last_success_frame== -1){
                        last_success_frame=i;
                    } else{
                        temp_waiting.push_back(condensed_state_frames[last_success_frame][3]);
                        temp_waiting.push_back(condensed_state_frames[i][1]);
                        temp_waiting.push_back(condensed_state_frames[i][1]-condensed_state_frames[last_success_frame][3]);
                        temp_waiting.push_back(condensed_state_frames[i+1][0]);
                        waiting_stat.push_back(temp_waiting);
                        temp_waiting.clear();
                        last_success_frame=i;
                    }
                }
            }
        } else{//Diffused out of boundary
            last_success_frame=-1;
        }
    }
    return waiting_stat;
}

int calc_distribution(std::vector<std::vector<int>> timelag_stat, std::vector<double> &time_distribution_array, int window){

    for (int i=0;i<timelag_stat.size();i++){
        if (timelag_stat[i][2]<max_transition_frames){
            time_distribution_array[timelag_stat[i][2]/window]++;
        }
    }
    return 0;
}
