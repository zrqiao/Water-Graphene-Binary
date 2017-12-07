//
// Created by utena on 17-10-27.
//
int Calc_Rotational_Relaxation_Correlation_Function(char* dist,int temperature,int start_nc,int end_nc) {
    CreatDir("Analysis/Rotational_Relaxation/OH_Bond");
    CreatDir("Analysis/Rotational_Relaxation/Dipole");
    CreatDir("Analysis/Rotational_Relaxation/Normal");
    typedef std::vector<double>::size_type index;
    std::cout << "program to calculate the Rotational Relaxation Correlation Function by clusters" << "\n" << std::endl;
    char out_name[64];
    amber_parm parm_name(name_parm7);
    std::vector <index0> O_WAT_id = parm_name.id_by_type("OW");
    int num[5];
    int WAT_ID_to_i[12000];
    for (int i = 0; i < 12000; i++) WAT_ID_to_i[i] = -1;
    for (int i = 0; i < O_WAT_id.size(); i++) WAT_ID_to_i[O_WAT_id[i]] = i;
    std::vector<int> Cluster_IDs[O_WAT_id.size()]={std::vector<int>()};
    std::ifstream cluster_infile;

}

std::vector<std::vector<int>>  Reduce(std::vector<std::vector<int>> frame_sta)
{
    typedef std::vector<double>::size_type index;
    std::vector<std::vector<int>> condensed_jump_frames;
    std::vector<int> temp_condensed_frames;//0-Cluster_size 1-startframe 2-framelength 3-endframe
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

std::vector<std::vector<int>>  pick_frame(int frame_start, int frame_end, std::vector<double>  XYZ_limit, int o_wat_id, int dens_ID) {
    std::vector<std::vector<int>> frame_in_graphene;
    std::vector<std::vector<int>> condensed_state_frames;
    char name_parm7[64];
    sprintf(name_parm7,"WAT_%d/density_dis9a5_WAT%d.parm7",WAT_NUM[dens_ID],WAT_NUM[dens_ID]);
    int frame_in_total_nc = 0;
    for(int nc = start_nc; nc  != end_nc+1; ++nc)
    {
        char name_nc[64];
        sprintf(name_nc, "WAT_%d/density_dis9a5_WAT%d_%d.nc", WAT_NUM[dens_ID],WAT_NUM[dens_ID],nc);
        amber_parm parm_name(name_parm7);
        nctraj nc_data(name_nc);
        std::vector<double>  o_coor;
        int frame_this_nc = nc_data.frames_number();
        for(int frame =0; frame < frame_this_nc; frame+=dt)
        {
            o_coor = nc_data.atom_coordinate(frame, o_wat_id);
            if( o_coor[2] < XYZ_limit[5] && o_coor[2] > XYZ_limit[4] && o_coor[1] < XYZ_limit[3] && o_coor[1] > XYZ_limit[2]
                && o_coor[0] < XYZ_limit[1] && o_coor[0] > XYZ_limit[0]) {//核心状态判断
                std::vector<int>  temp;
                temp.push_back(frame_in_total_nc+frame);
                if(o_coor[2] > ((XYZ_limit[5]*80+XYZ_limit[4]*80)/160 + width_transition)) {
                    temp.push_back(1);
                }

                frame_in_graphene.push_back(temp);
            }
        }
        frame_in_total_nc = frame_in_total_nc + frame_this_nc;
    }

    condensed_state_frames=Reduce(frame_in_graphene);
    frame_in_graphene.clear();
    return condensed_state_frames;
}

std::vector<std::vector<int>>  calc_exchange( std::vector<std::vector<int>> condensed_state_frames){
    std::vector<std::vector<int>> exchange_stat;
    std::vector<int> temp_jump;// 0-startframe 1-endframe 2-transitiontime 3-direction
    for (int i=1;i<condensed_state_frames.size()-1;i++){
        if (condensed_state_frames[i][0]==0){
            if(condensed_state_frames[i-1][3]+dt==condensed_state_frames[i][1]&&condensed_state_frames[i][3]+dt==condensed_state_frames[i+1][1]){//修改
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