//
// Created by utena on 17-10-27.
//



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
                std::vector<int>  bbb;
                bbb.push_back(frame_in_total_nc+frame);
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
        }
        frame_in_total_nc = frame_in_total_nc + frame_this_nc;
    }
    //~ if(frame_in_graphene.size() != 0 )
    //~ {
    //~ while(frame_in_graphene[0][1] == 55)
    //~ {
    //~ frame_in_graphene.erase(frame_in_graphene.begin());
    //~ }
    //~ }
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