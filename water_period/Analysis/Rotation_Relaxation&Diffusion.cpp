//
// Created by utena on 17-10-27.
//
#include "Environment& Algorithm.hpp"
#include <string>
#define max_transition_time 20000
#define tot_frame 20000
#define frame_segment 50
#define coarse_grain_cutoff 8
#define max_cluster_size 24
using std::cout;
using std::vector;
using std::ifstream;
using std::cin;
using std::string;
using std::endl;
using std::setw;

int  Reduce(std::vector<std::vector<int>> frame_sta, std::vector<std::vector<int>> &condensed_frames);
int pick_frame(int frame_zeropoint, std::vector<int> Cluster_Size_by_WATID_frame,std::vector<std::vector<int>>& condensed_state_frames);

int Calc_Rotational_Relaxation_Correlation_Function(char* dist,int temperature,int start_nc,int end_nc) {
    CreatDir("Analysis/Rotational_Relaxation/OH_Bond");
    CreatDir("Analysis/Rotational_Relaxation/Dipole");
    CreatDir("Analysis/Rotational_Relaxation/Normal");
    CreatDir("Analysis/Cluster_size");
    std::cout << "program to calculate the Rotational Relaxation Correlation Function by clusters" << "\n" << std::endl;
    cout<<"Loading trajectories..."<<endl;
    amber_parm parm_name(name_parm7);
    std::vector <index0> O_WAT_id = parm_name.id_by_type("OW");
    int tot_wat_num=O_WAT_id.size();
    std::vector<int> Cluster_IDs[tot_wat_num]={std::vector<int>()};
    std::vector<int> Cluster_Size_by_WATID_frame[tot_wat_num]={std::vector<int>()};
    std::vector<double> OH1_vec(3,0),OH2_vec(3,0),Dip_vec(3,0),Nor_vec(3,0);
    std::vector<double> OH1_vec_0(3,0),OH2_vec_0(3,0),Dip_vec_0(3,0),Nor_vec_0(3,0);
    std::vector<double> O_coor, O1_coor,O2_coor, H1O_coor, H2O_coor,H1_coor, H2_coor,O_coor_0;
    double RCF_sum_OH [max_cluster_size/coarse_grain_cutoff][frame_segment]={0};
    double RCF_sum_Dip[max_cluster_size/coarse_grain_cutoff][frame_segment]={0};
    double RCF_sum_Nor[max_cluster_size/coarse_grain_cutoff][frame_segment]={0};
    double Diffusion_sum[max_cluster_size/coarse_grain_cutoff][frame_segment]={0};
    int tot_count[max_cluster_size/coarse_grain_cutoff][frame_segment]={0};
    int buff[5];
    int WAT_ID_to_i[12000];
    for (int i = 0; i < 12000; i++) WAT_ID_to_i[i] = -1;
    for (int i = 0; i < O_WAT_id.size(); i++) WAT_ID_to_i[O_WAT_id[i]] = i;
    ifstream cluster_infile;
    char infile_name[64];
    char name_nc[64];
    for (int nc=start_nc;nc<=end_nc;nc++){
        sprintf(infile_name, "Analysis/Cluster_Data/nc_%02d.dat", nc);
        cluster_infile.open(infile_name);
        //sprintf(name_nc, "density_dis%s_WAT1150_NVT_%d_unwrapped_%d.nc",dist,temperature,nc);
        //nctraj nc_data(name_nc);
        //int tot_frame=nc_data.frames_number();
        string mark;
        int frame=(nc-start_nc)*tot_frame;//爲易於迭代以開始nc作爲0時間 注意輸出時調整時間零點
        while (getline(cluster_infile,mark)){
            if (frame%1000==0) cout<<mark<<endl;
            for (int WAT_ID_ref=0;WAT_ID_ref<tot_wat_num;WAT_ID_ref++){
                Cluster_IDs[WAT_ID_ref].clear();
            }
            //cout<<frame<<endl;
            for (int WAT_ID_ref=0;WAT_ID_ref<tot_wat_num;WAT_ID_ref++){
                cluster_infile>>buff[0]>>buff[1];
                int WAT_ID=buff[0],Cluster_ID=buff[1];
                Cluster_IDs[Cluster_ID].push_back(WAT_ID_to_i[WAT_ID]);
                //cout<<WAT_ID<<endl;
            }
            for (int WAT_ID_ref=0;WAT_ID_ref<O_WAT_id.size();WAT_ID_ref++){
                Cluster_Size_by_WATID_frame[WAT_ID_ref].push_back(Cluster_IDs[WAT_ID_ref].size());
            }
            frame++;
            getline(cluster_infile,mark);
        }
        cluster_infile.close();
    }
    /*for (int WAT_ID_ref = 0; WAT_ID_ref != tot_wat_num; ++WAT_ID_ref) {
        char outfile_name[64];
        std::ofstream tracing_out;
        sprintf(outfile_name, "Analysis/Cluster_size/WAT_%d_nc_%02d-%02d.dat",O_WAT_id[WAT_ID_ref], start_nc,end_nc);
        tracing_out.open(outfile_name);
        for (int t=0;t!=Cluster_Size_by_WATID_frame[0].size();t++){
            tracing_out<<t*0.25<<setw(10)<<Cluster_Size_by_WATID_frame[WAT_ID_ref][t]<<endl;
        }
        tracing_out.close();
    }*/
    cout<<"Loading complete"<<endl;
    int nc_id=start_nc,frame_id=-1,frame_r=-1;
    sprintf(name_nc, "density_dis%s_WAT1150_NVT_%d_unwrapped_%d.nc", dist, temperature, nc_id);
    nctraj nc_data(name_nc);
    for (int WAT_ID_ref = 0; WAT_ID_ref != tot_wat_num; ++WAT_ID_ref) {//calculate RCF
        cout << "water calculate and id: " << setw(8) << O_WAT_id[WAT_ID_ref] << endl;
        std::vector<std::vector<int>> pickframe;
        pick_frame(start_nc * tot_frame, Cluster_Size_by_WATID_frame[WAT_ID_ref],pickframe);
        //for (index state_id = 0; state_id != pickframe.size(); ++state_id) {
          //  std::cout << pickframe[state_id][0] << "   " << pickframe[state_id][1] << "   " << pickframe[state_id][2]
            //          << "   "
              //        << pickframe[j][state_id] << std::endl;
        //}
        int Target_ID = O_WAT_id[WAT_ID_ref];
        if (pickframe.size() != 0) {
            for (int j = 0; j < pickframe.size(); j++) {
                //cout<<"start frame:"<<pickframe[j][1]<<" size:"<<pickframe[j][0]<<endl;
                for (int frame_raw = pickframe[j][1]; frame_raw < pickframe[j][3]; frame_raw++) {
                    if (nc_id!=floor(frame_raw / tot_frame)) {
                        nc_id = floor(frame_raw / tot_frame);
                        sprintf(name_nc, "density_dis%s_WAT1150_NVT_%d_unwrapped_%d.nc", dist, temperature, nc_id);
                        nctraj nc_data(name_nc);
                    }
                    frame_id = frame_raw - nc_id * tot_frame;
                    O_coor = nc_data.atom_coordinate(frame_id, Target_ID);
                    H1_coor = nc_data.atom_coordinate(frame_id, Target_ID + 1);
                    H2_coor = nc_data.atom_coordinate(frame_id, Target_ID + 2);
                    OH1_vec = vector_calc::vector_minus(H1_coor, O_coor);
                    OH2_vec = vector_calc::vector_minus(H2_coor, O_coor);
                    Dip_vec = vector_calc::vector_minus(
                            vector_calc::vector_divide(vector_calc::vector_add(H1_coor, H2_coor), 2), O_coor);
                    Nor_vec = vector_calc::cross_product(OH2_vec, OH1_vec);
                    if ((frame_raw - pickframe[j][1]) % frame_segment == 0) {//时间零点
                        OH1_vec_0 = OH1_vec;
                        OH2_vec_0 = OH2_vec;
                        Dip_vec_0 = Dip_vec;
                        Nor_vec_0 = Nor_vec;
                        O_coor_0 = O_coor;
                    }
                    frame_r = (frame_raw - pickframe[j][1]) % frame_segment;
                    //std::cout<<second_order_legendre_function(vector_calc::vector_angle(OH1_vec,OH1_vec_0))<<std::endl;
                    tot_count[pickframe[j][0]][frame_r] += 1;
                    RCF_sum_Dip[pickframe[j][0]][frame_r] += factor_cos(vector_calc::vector_angle(Dip_vec, Dip_vec_0));
                    RCF_sum_Nor[pickframe[j][0]][frame_r] += factor_cos(vector_calc::vector_angle(Nor_vec, Nor_vec_0));
                    RCF_sum_OH[pickframe[j][0]][frame_r] += factor_cos(vector_calc::vector_angle(OH1_vec, OH1_vec_0));
                    RCF_sum_OH[pickframe[j][0]][frame_r] += factor_cos(vector_calc::vector_angle(OH2_vec, OH2_vec_0));
                    Diffusion_sum[pickframe[j][0]][frame_r] += pow(
                            vector_calc::vector_module(vector_calc::vector_minus(O_coor, O_coor_0)), 2);
                }
            }
            //calc_distribution(jump_stat, transitionpath_time_distribution);
            pickframe.clear();
            //vector<vector<int>>().swap(pickframe);
            int j;
            char out_name[64];
            std::ofstream distribution_out;
            sprintf(out_name, "Analysis/Rotational_Relaxation/OH_Bond/nc%02d-%02d.dat", start_nc,end_nc);
            distribution_out.open(out_name);
            for (int i = 0; i != frame_segment; i++) {
                distribution_out << i * 0.25 << std::setw(15);
                for (int j=0;j!=floor(max_cluster_size/coarse_grain_cutoff);j++){
                    distribution_out << RCF_sum_OH[j][i] / (tot_count[j][i] * 2) << std::setw(15);
                }
                distribution_out << std::endl;
            }
            distribution_out.close();
            sprintf(out_name, "Analysis/Rotational_Relaxation/Dipole/nc%02d-%02d.dat", start_nc,end_nc);
            distribution_out.open(out_name);
            for (int i = 0; i != frame_segment; i++) {
                distribution_out << i * 0.25 << std::setw(15);
                for (j=0;j!=floor(max_cluster_size/coarse_grain_cutoff);j++){
                    distribution_out << RCF_sum_Dip[j][i] / (tot_count[j][i]) << std::setw(15);
                }
                distribution_out << std::endl;
            }
            distribution_out.close();
            sprintf(out_name, "Analysis/Rotational_Relaxation/Normal/nc%02d-%02d.dat", start_nc,end_nc);
            distribution_out.open(out_name);
            for (int i = 0; i != frame_segment; i++) {
                distribution_out << i * 0.25 << std::setw(15);
                for (j=0;j!=floor(max_cluster_size/coarse_grain_cutoff);j++){
                    distribution_out << RCF_sum_Nor[j][i] / (tot_count[j][i]) << std::setw(15);
                }
                distribution_out << std::endl;
            }
            distribution_out.close();
            sprintf(out_name, "Analysis/Diffusion/nc%02d-%02d.dat", start_nc,end_nc);
            distribution_out.open(out_name);
            for (int i = 0; i != frame_segment; i++) {
                distribution_out << i * 0.25 << std::setw(15);
                for (j=0;j!=floor(max_cluster_size/coarse_grain_cutoff);j++){
                    distribution_out << Diffusion_sum[j][i] / (tot_count[j][i]) << std::setw(15);
                }
                distribution_out << std::endl;
            }
            distribution_out.close();
        }


    }

    return 0;
}

int  Reduce(std::vector<std::vector<int>> frame_sta, std::vector<std::vector<int>> &condensed_frames)
{
    //typedef std::vector<double>::size_type index;
    std::vector<int> temp_condensed_frames;//0-Cluster_size 1-startframe 2-framelength 3-endframe
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

int pick_frame(int frame_zeropoint, std::vector<int> Cluster_Size_by_WATID_frame,std::vector<std::vector<int>>& condensed_state_frames) {
    std::vector<std::vector<int>> frame_in_graphene;
    for(int frame =0; frame < Cluster_Size_by_WATID_frame.size()-4; frame+=dt)
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
    Reduce(frame_in_graphene,condensed_state_frames);
    frame_in_graphene.clear();
    return 0;
}

std::vector<std::vector<int>>  calc_exchange( std::vector<std::vector<int>> condensed_state_frames){
    std::vector<std::vector<int>> exchange_stat;
    std::vector<int> temp_jump;// 0-startframe 1-endframe 2-transitiontime 3-direction
    for (int i=1;i<condensed_state_frames.size()-1;i++){
        if (condensed_state_frames[i][0]==0){
            if(condensed_state_frames[i-1][3]+dt==condensed_state_frames[i][1]&&condensed_state_frames[i][3]+dt==condensed_state_frames[i+1][1]){//修改
                if (condensed_state_frames[i-1][0]!=condensed_state_frames[i+1][0]){
                    temp_jump.push_back(condensed_state_frames[i][1]);
                    temp_jump.push_back(condensed_state_frames[i][3]);
                    temp_jump.push_back(condensed_state_frames[i][2]);
                    temp_jump.push_back(condensed_state_frames[i+1][0]);
                    exchange_stat.push_back(temp_jump);
                    temp_jump.clear();
                }
            }
        }
    }
    return exchange_stat;
}
int calc_distribution(std::vector<std::vector<int>> jump_stat,std::vector<double> &transitiontime_distribution){

    for (int i=0;i<jump_stat.size();i++){
        if (jump_stat[i][2]<max_transition_time){
            transitiontime_distribution[jump_stat[i][2]/dt]++;
        }
    }
    return 0;
}