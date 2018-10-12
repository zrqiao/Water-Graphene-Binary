#include "amber_netcdf.hpp"
#include "amber_parm_1.0.hpp"
#include "vector_calc.h"
//~ #include <string.h>
//~ #include <fstream>
//~ #include <iostream>
//~ #include <stdlib.h>
//~ #include <stdio.h>


#define z_axis_modify  2
#define deviation_y_center 20.0
#define deviation_x_center 20.0
#define start_nc 16
#define end_nc  16
#define z_step 65
#define angle_number 90
#define name_parm7 "density_dis9a5.parm7" 
//~ #define name_nc "water_ion_graphene_10a5"

int main()
{
	typedef std::vector<double>::size_type index;
	std::cout << "program to calculate the z axis distribution of dangling"<< "\n"<< std::endl;
	
	int dangling_distribution[z_step][angle_number]={0};
	
	for(int nc = start_nc; nc  != end_nc+1; ++nc)
	{
        char name_nc[64];
        sprintf(name_nc, "density_dis9a5_%d.nc",nc);
		amber_parm parm_name(name_parm7);
        nctraj nc_data(name_nc);
        printf("nc_file = %d\n",nc);
        int total_frame = nc_data.frames_number();
		std::cout << "total frame: " << total_frame << std::endl;
	
	
	    std::vector<index> O_WAT_id = parm_name.id_by_type("OW");	
     
       for(int frame =0; frame != total_frame; ++frame)
       //~ for(int frame =1000; frame != 1001; ++frame)
       {
		   if((frame %100) == 0)
		   {
		   std::cout << "frame_now :            " << frame <<std::endl;
	       }
	     
	       double C_z_coor_sum_1 = 0, C_z_coor_sum_2 = 0;
           double C_z_coor_average_1, C_z_coor_average_2;	
	       double C_x_coor_sum = 0, C_y_coor_sum = 0;
           double C_x_coor_center, C_y_coor_center;        
           double Z_UP,Z_DOWN,Y_UP,Y_DOWN,X_UP,X_DOWN;	      
            
	        for( index C_index = 0 ; C_index != 1400 ; ++C_index)
	        {
		           C_z_coor_sum_1 += nc_data.atom_coordinate(frame, C_index)[2];
		           C_z_coor_sum_2 += nc_data.atom_coordinate(frame, (1400+C_index))[2];
		           C_y_coor_sum += nc_data.atom_coordinate(frame, C_index)[1];
   		           C_x_coor_sum += nc_data.atom_coordinate(frame, C_index)[0];		
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
	          
	         std::vector<index> O_WAT_IN_C_id;       
	         std::vector<double> O_coor,O_coor_initial;
	
		      for(index i = 0; i != O_WAT_id.size(); ++i)
			  {
			      O_coor = nc_data.atom_coordinate(frame, O_WAT_id[i]);
				   if( O_coor[2] < Z_UP && O_coor[2] > Z_DOWN && O_coor[0] < X_UP && O_coor[0] > X_DOWN && O_coor[1] < Y_UP && O_coor[1] > Y_DOWN)
				   {
						 O_WAT_IN_C_id.push_back(O_WAT_id[i]);
					}			
		       }

		      for(index i =0; i != O_WAT_IN_C_id.size(); ++i)
		      {	  
						  //~ double distance;
						  std::vector<double> orientation(3,0);
						  std::vector<double> z_axis={0,0,1};
						  orientation[0]=nc_data.atom_coordinate(frame,O_WAT_IN_C_id[i]+1)[0]-nc_data.atom_coordinate(frame,O_WAT_IN_C_id[i])[0];
						  orientation[1]=nc_data.atom_coordinate(frame,O_WAT_IN_C_id[i]+1)[1]-nc_data.atom_coordinate(frame,O_WAT_IN_C_id[i])[1];	
						  orientation[2]=nc_data.atom_coordinate(frame,O_WAT_IN_C_id[i]+1)[2]-nc_data.atom_coordinate(frame,O_WAT_IN_C_id[i])[2];						  					  
						  //~ distance   = sqrt(pow((nc_data.atom_coordinate(frame,O_WAT_IN_C_id[i])[0]- nc_data.atom_coordinate(frame,O_WAT_IN_C_id[i]+1)[0]),2) +
						                         //~ pow((nc_data.atom_coordinate(frame,O_WAT_IN_C_id[i])[1]- nc_data.atom_coordinate(frame,O_WAT_IN_C_id[i]+1)[1]),2) +
						                         //~ pow((nc_data.atom_coordinate(frame,O_WAT_IN_C_id[i])[2]- nc_data.atom_coordinate(frame,O_WAT_IN_C_id[i]+1)[2]),2));	

							  double angle;
							  angle = vector_calc::vector_angle(z_axis,orientation);
							  //~ std::cout<< angle<<std::endl;
							  //~ std::cout<<(Z_UP-Z_DOWN)/z_step<<std::endl;							  
							  //~ std::cout<<floor((nc_data.atom_coordinate(frame,O_WAT_IN_C_id[i]+1)[2]-Z_DOWN)/((Z_UP-Z_DOWN)/z_step))<<std::endl;
							   ++dangling_distribution[int(floor((nc_data.atom_coordinate(frame,O_WAT_IN_C_id[i]+1)[2]-Z_DOWN)/((Z_UP-Z_DOWN)/z_step)))][int(floor((180.0*acos(angle)/vector_calc::PI)/(180.0/angle_number)))];
							  					      				  
			  }	 
	     }       // frame end     
	}      // nc end

    std::ofstream outfile;
    outfile.open("dangling_distribution_unnormal.dat");
    
    for(int i = 0; i!=z_step; ++i)
    //~ for(int i = 0; i!=1; ++i) 
    {
		for(int j =0; j != angle_number; ++j)
		{
			std::cout << sin((j+0.5)*(vector_calc::PI/angle_number))<<std::endl;
			//~ std::cout << sin(acos((j+0.5)*(180.0/angle_number)))<<std::endl;
			//~ outfile<<std::setw(10)<< dangling_distribution[i][j]/sin((j+0.5)*(vector_calc::PI/angle_number));
	        outfile<<std::setw(10)<< dangling_distribution[i][j];		
		}
		outfile<<std::endl;
	}

	return 0;
}


