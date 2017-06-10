#include "amber_netcdf.hpp"
#include "amber_parm_1.0.hpp"
#include "vector_calc.h"


//~ #define z_axis_modify  2
#define width_transition 1.3
#define deviation_y_center 26
#define deviation_x_center 26
#define start_nc 12
#define end_nc  16
//~ #define relaxation_time 25
#define name_parm7 "density_dis9a5.parm7" 
#define jump_time 5000
//~ #define name_nc "water_ion_graphene_10a5"

std::vector<std::vector<int>>  pick_frame(int frame_start, int  frame_end, std::vector<double>   XYZ_limit, int o_wat_id);
std::vector<int>  jump_count(std::vector<std::vector<int>> frame_sta,std::vector<int> count_statistics);

int main()
{
	typedef std::vector<double>::size_type index;
	std::cout << "program to calculate the jump time distribution"<< "\n"<< std::endl;
	std::cout << "calculate XYZ_limit" << std::endl;
	
    double C_z_coor_sum_1 = 0, C_z_coor_sum_2 = 0;
    double C_z_coor_average_1, C_z_coor_average_2;	
	double C_x_coor_sum = 0, C_y_coor_sum = 0;
    double C_x_coor_center, C_y_coor_center;     
    double Z_UP,Z_DOWN,Y_UP,Y_DOWN,X_UP,X_DOWN;
	
	amber_parm parm_nam(name_parm7);
	nctraj data_nc("density_dis9a5_15.nc");

	for( index C_index = 0 ; C_index != 1400 ; ++C_index)
	{
		       C_z_coor_sum_1 += data_nc.atom_coordinate(0, C_index)[2];
		       C_z_coor_sum_2 += data_nc.atom_coordinate(0, (1400+C_index))[2];
		       C_y_coor_sum += data_nc.atom_coordinate(0, C_index)[1];
   		       C_x_coor_sum += data_nc.atom_coordinate(0, C_index)[0];		
	 }		    
	 C_z_coor_average_1 = C_z_coor_sum_1/1400;
	 C_z_coor_average_2 = C_z_coor_sum_2/1400;
	 C_y_coor_center = C_y_coor_sum/1400;
	 C_x_coor_center = C_x_coor_sum/1400;
	            
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
	
	 
    std::vector<std::vector<int>> pickframe;		
    std::vector<index> O_WAT_id = parm_nam.id_by_type("OW");	
	std::cout<< "total water:  " << O_WAT_id.size() <<"\n" <<std::endl;
	
    std::vector<int> count_all_water(jump_time,0);
	std::vector<int> wat_jump_count;
	std::vector<std::vector<int>> wat_jump_distribution;
	
	
	for(index i = 0; i != O_WAT_id.size(); ++i)
	//~ for(index i = 0; i != 1; ++i)	
    {
		std::cout << "water calculate and id: " << std::setw(8)<< i+1 <<std::setw(8)<<O_WAT_id[i]<< std::endl;
		pickframe = pick_frame(start_nc,end_nc,xyz,O_WAT_id[i]);
		//~ std::cout<<"frame size: " << pickframe.size()<<std::endl;
		for(index j =0; j != pickframe.size(); ++j) std::cout << pickframe[j][0] <<"   "<<pickframe[j][1] << std::endl;
		if(pickframe.size() != 0)
		{
		    count_all_water = jump_count(pickframe,count_all_water);
		    
		}
    }
	
	//~ std::ofstream outfile;

	
	//~ outfile.open("count_jump_time_with_transition");

	for(index i = 2; i < count_all_water.size(); i += 4)
	{
		//~ outfile <<std::setw(15) <<i<<std::setw(15)<<count_all_water[i] << std::endl;
		std::cout <<std::setw(15) <<i+1.5<<std::setw(15)<<(count_all_water[i]+count_all_water[i+1]+count_all_water[i+2]+count_all_water[i+3]) << std::endl;
     }
     //~ outfile.close();      
	//~ return 0;
}



std::vector<std::vector<int>>  pick_frame(int frame_start, int frame_end, std::vector<double>  XYZ_limit, int o_wat_id)
{
	std::vector<std::vector<int>> frame_in_graphene;
	int frame_in_total_nc = 0;
    for(int nc = start_nc; nc  != end_nc+1; ++nc)
    {
        char name_nc[64];
        sprintf(name_nc, "density_dis9a5_%d.nc",nc);
		amber_parm parm_name(name_parm7);
        nctraj nc_data(name_nc);
        std::vector<double>  o_coor;
        int frame_this_nc = nc_data.frames_number();
        for(int frame =0; frame != frame_this_nc; ++frame)
        {
			o_coor = nc_data.atom_coordinate(frame, o_wat_id);
			if( o_coor[2] < XYZ_limit[5] && o_coor[2] > XYZ_limit[4] && o_coor[1] < XYZ_limit[3] && o_coor[1] > XYZ_limit[2] && o_coor[0] < XYZ_limit[1] && o_coor[0] > XYZ_limit[0]) 
			{
				std::vector<int>  bbb;
				bbb.push_back(frame_in_total_nc+frame);
				if(o_coor[2] > ((XYZ_limit[5]+XYZ_limit[4])/2 + width_transition))
				{
					bbb.push_back(1);
				 }
				 else if(o_coor[2] < ((XYZ_limit[5]+XYZ_limit[4])/2 - width_transition))
				 {
					 bbb.push_back(0);
			     }
			     else
			     {
					 bbb.push_back(55);
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
	return frame_in_graphene;
}	

std::vector<int>  jump_count(std::vector<std::vector<int>> frame_sta,std::vector<int> count_statistics)
{
	typedef std::vector<double>::size_type index;
	int start1 = frame_sta[0][1];
	int start2 ;
	int start2_logic=0;
	int count=1;																																			//~ std::cout <<count<<std::endl;
    for(index i=0; i != (frame_sta.size()-1); ++i)
	{
		if(frame_sta[i+1][0] ==  (frame_sta[i][0]+1))
		{
			if(start1 == 55)
			{
				start1= frame_sta[i+1][1];
			} 
			if(start1 != 55)
			{
			    if(start2_logic == 0)
			    {
			        if((frame_sta[i+1][1] != start1) && (frame_sta[i+1][1] != 55))
			        {
				        start2  =  frame_sta[i+1][1];
				        start2_logic = 1;
				        count = 0;
				        //~ std::cout << "frame: "<<frame_sta[i+1][0]<< " start2_logic: " << start2_logic <<std::endl;
			        }
			    }
			}
			if(frame_sta[i+1][1] == frame_sta[i][1])
			{
				count += 1;			
				    //~ std::cout << "frame: "<<frame_sta[i+1][0]<< "  count+1 for frame_sta[i+1][1] == frame_sta[i][1]"<<std::endl;
			}
			else if(frame_sta[i+1][1] ==  55)
			{
				count += 1;
				    //~ std::cout << "frame: "<<frame_sta[i+1][0]<< "  count+1 for frame_sta[i+1][1] == 55"<<std::endl;																														//~ std::cout <<count<<std::endl;
			}	
			else if((start2_logic == 0) && (frame_sta[i+1][1] == start1))
			{ 
			    count += 1;
				    //~ std::cout << "frame: "<<frame_sta[i+1][0]<< "  count+1 for (start2_logic == 0) && (frame_sta[i+1][1] == start1)"<<std::endl;		
			    
			}																		
			else if((start2_logic == 1) && (frame_sta[i+1][1] == start2))
			{
				count += 1;			
				    //~ std::cout << "frame: "<<frame_sta[i+1][0]<< "  count+1 for (start2_logic == 1) && (frame_sta[i+1][1] == start2)"<<std::endl;																																											//~ std::cout <<count<<std::endl;
			}
			else
			{	
				//~ std::cout<<"statistics: "<<count<<"  next start: "<< frame_sta[i+1][0]<<std::endl;
				count += 1;
				count_statistics[count] += 1 ;
				start2 = frame_sta[i+1][1];
				count = 1;
				std::cout << "frame: "<<frame_sta[i+1][0]<< "  count=1 for start1 = 1"<<std::endl;	
			}
	    }
	    else if(frame_sta[i+1][0] !=  (frame_sta[i][0]+1))
	    {
			//~ ++i;
			//~ start = frame_sta[i+1][1];
			//~ while(((i+2) < (frame_sta.size()-2))   && (frame_sta[i+2][0] == (frame_sta[i+1][0]+1))) ++i;
			//~ while(((i+1) <  (frame_sta.size()-1))  &&  (frame_sta[i+1][1] == 55) )
			//~ {
				//~ ++i;
			//~ }
			start1 = frame_sta[i+1][1];
			start2_logic = 0;
			count = 1;
				//~ std::cout << "frame: "<<frame_sta[i+1][0]<< "  count=1 for frame_sta[i+1][0] !=  (frame_sta[i][0]+1)"<<std::endl;																																																												//~ std::cout <<count<<std::endl;
		}
		
    }
    
    return count_statistics;
}

