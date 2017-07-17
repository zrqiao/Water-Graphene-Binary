//
// Created by utena on 17-7-16.
//

#include "amber_netcdf.hpp"
#include "amber_parm_1.0.hpp"
#include "vector_calc.h"
//~ #define z_axis_modify  2
#define width_transition 0.45
#define deviation_y_center 26
#define deviation_x_center 26
#define start_nc 37
#define end_nc  100
//~ #define relaxation_time 25
#define name_parm7 "density_dis9a5.parm7"
#define jump_time 5000
#define max_transition_time 2000
#define dt 4
typedef std::vector<double>::size_type index;
//~ #define name_nc "water_ion_graphene_10a5"
std::vector<std::vector<int>>  pick_frame(int frame_start, int  frame_end, std::vector<double>   XYZ_limit, int o_wat_id);
std::vector<std::vector<int>>  calc_jump(std::vector<std::vector<int>> condensed_frames);
int calc_distribution(std::vector<std::vector<int>> jump_stat,std::vector<double> &transitiontime_distribution);
int main()
{
    std::vector<int*> jump_coor;
    std::vector<int> recrossing_id;
    std::vector<std::vector <int>> recrossing_stat;
    std::cout << "program to calculate recrossing" << "\n" << std::endl;
    std::vector<int> jump_time_distribution(jump_time,0);
    std::cout << "calculate XYZ_limit" << std::endl;
    std::vector<double> recrossing_time_distribution(max_transition_time,0);
    std::vector<double> recrossing_count_distribution(50,0);
    std::vector<std::vector<int>> jump_stat;

    double C_z_coor_sum_1 = 0, C_z_coor_sum_2 = 0;
    double C_z_coor_average_1, C_z_coor_average_2;
    double C_x_coor_sum = 0, C_y_coor_sum = 0;
    double C_x_coor_center, C_y_coor_center;
    double Z_UP,Z_DOWN,Y_UP,Y_DOWN,X_UP,X_DOWN;
    for (int nc=start_nc;nc<=end_nc;nc++) {
        amber_parm parm_nam(name_parm7);
        char name_nc[64];
        sprintf(name_nc, "nc/density_dis9a5_%d.nc",nc);
        nctraj data_nc(name_nc);
        for (index C_index = 0; C_index != 1400; ++C_index) {
            C_z_coor_sum_1 += data_nc.atom_coordinate(0, C_index)[2];
            C_z_coor_sum_2 += data_nc.atom_coordinate(0, (1400 + C_index))[2];
            C_y_coor_sum += data_nc.atom_coordinate(0, C_index)[1];
            C_x_coor_sum += data_nc.atom_coordinate(0, C_index)[0];
        }
    }
    C_z_coor_average_1 = C_z_coor_sum_1/(1400*(end_nc-start_nc+1));
    C_z_coor_average_2 = C_z_coor_sum_2/(1400*(end_nc-start_nc+1));
    C_y_coor_center = C_y_coor_sum/(1400*(end_nc-start_nc+1));
    C_x_coor_center = C_x_coor_sum/(1400*(end_nc-start_nc+1));

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

    std::vector<double> xyz;
    xyz.push_back(X_DOWN);
    xyz.push_back(X_UP);
    xyz.push_back(Y_DOWN);
    xyz.push_back(Y_UP);
    xyz.push_back(Z_DOWN);
    xyz.push_back(Z_UP);

    amber_parm parm_nam(name_parm7);
    std::vector<std::vector<int>> pickframe;
    std::vector<index> Target_ID = parm_nam.id_by_type("OW");
    std::cout<< "total water:  " << Target_ID.size() <<"\n" <<std::endl;

    std::vector<int> wat_jump_count;
    std::vector<std::vector<int>> wat_jump_distribution;
    std::ifstream infile;
    infile.open("jump/cutoff_1a2/recrossing_index_start_finish_down_2");
    double num[4];
    while (!infile.fail()) {
        infile >> num[0] >> num[1]>>num[2]>>num[3];
        std::cout<<num[0]<<std::setw(10)<<num[1]<<std::endl;
        int frame_r_st = num[1]-5;
        int frame_r_ed = num[2]+5;
        int Target_ID = num[0];
        int Direction = num[3];
        pickframe = pick_frame(frame_r_st, frame_r_ed, xyz, Target_ID);
        //~ std::cout<<"frame size: " << pickframe.size()<<std::endl;
        /*for (index j = 0; j != pickframe.size(); ++j)
            std::cout << pickframe[j][0] << "   " << pickframe[j][1] << "   " << pickframe[j][2] << "   "
                      << pickframe[j][3] << std::endl;*/
        if (pickframe.size() != 0) {
            jump_stat = calc_jump(pickframe);

            std::cout<<jump_stat.size()<<std::endl;
            if (jump_stat.size()==1) {
                for (int j = 0; j < jump_stat.size(); j++) {
                    recrossing_stat.push_back(jump_stat[j]);
                    recrossing_id.push_back(Target_ID);
                }
            }
            recrossing_count_distribution[jump_stat.size()]++;
            calc_distribution(jump_stat, recrossing_time_distribution);
            jump_stat.clear();
            pickframe.clear();
        }
    }
    std::ofstream outfile2;
    outfile2.open("jump/cutoff_1a2/recrossing_index_start_finish_down_3");
    for (int j = 0; j < recrossing_id.size(); j++) {
        outfile2 << recrossing_id[j] << std::setw(10) << recrossing_stat[j][0] << std::setw(10) << recrossing_stat[j][1]
                 << std::setw(10) << recrossing_stat[j][3] << std::endl;
    }
    outfile2.close();
    std::ofstream outfile;
    std::ofstream outfile1;
    outfile1.open("jump/cutoff_1a2/recrossing_count_in_jump_distribution_3");
    outfile.open("jump/cutoff_1a2/recrossing_time_distribution_3");
    for (index i = 0; i < recrossing_time_distribution.size(); i += 1) {
        //~ outfile <<std::setw(15) <<i<<std::setw(15)<<count_all_water[i] << std::endl;
        outfile << i * dt * 4 << std::setw(15) << recrossing_time_distribution[i] << std::setw(15) << std::endl;
    }
    outfile.close();
    for (index i = 0; i < 50; i += 1) {
        //~ outfile <<std::setw(15) <<i<<std::setw(15)<<count_all_water[i] << std::endl;
        outfile1 << i  << std::setw(15) << recrossing_count_distribution[i] << std::setw(15) << std::endl;
    }
    outfile1.close();
    return 0;
}


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
        if (i==frame_sta.size()-1){
            temp_condensed_frames.push_back(frame_sta[i][0]);
            condensed_jump_frames.push_back(temp_condensed_frames);
            temp_condensed_frames.clear();
        }
    }

    return condensed_jump_frames;
}


std::vector<std::vector<int>>  pick_frame(int frame_start, int frame_end, std::vector<double>  XYZ_limit, int o_wat_id) {
    std::vector<std::vector<int>> frame_in_graphene;
    std::vector<std::vector<int>> condensed_state_frames;
    std::vector<double> o_coor;
    char name_nc[64];
    int nc = start_nc;
    int total_frame = (start_nc) * 10000;
    sprintf(name_nc, "nc/density_dis9a5_%d.nc", nc);
    nctraj nc_data(name_nc);
    int totframe_nc=nc_data.frames_number();
    for(int frame_r =frame_start; frame_r <= frame_end; frame_r+=dt)
    {
        if (frame_r==frame_start || frame_r>=(total_frame+totframe_nc)){
            while (total_frame+totframe_nc <= frame_r) {
                total_frame += totframe_nc;
                nc++;
                sprintf(name_nc, "nc/density_dis9a5_%d.nc", nc);
                nctraj nc_data(name_nc);
                totframe_nc=nc_data.frames_number();
            }
        }
        sprintf(name_nc, "nc/density_dis9a5_%d.nc", nc);
        nctraj nc_data(name_nc);
        index frame = frame_r - total_frame ;
        o_coor = nc_data.atom_coordinate(frame, o_wat_id);

        std::vector<int>  bbb;
        bbb.push_back(frame_r);
        if(o_coor[2] > ((XYZ_limit[5]*80+XYZ_limit[4]*80)/160 + width_transition))
        {
            bbb.push_back(1);
        }
        else if(o_coor[2] < ((XYZ_limit[5]*80+XYZ_limit[4]*80)/160 - width_transition))
        {
            bbb.push_back(-1);
        }
        else
        {
            bbb.push_back(0);
        }
        frame_in_graphene.push_back(bbb);
    }
    condensed_state_frames=Reduce(frame_in_graphene);
    frame_in_graphene.clear();
    return condensed_state_frames;
}

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
int calc_distribution(std::vector<std::vector<int>> jump_stat,std::vector<double> &transitiontime_distribution){

    for (int i=0;i<jump_stat.size();i++){
        if (jump_stat[i][2]<max_transition_time){
            transitiontime_distribution[jump_stat[i][2]/dt]++;
        }
    }
    return 0;
}