#include "amber_netcdf.hpp"
#include "amber_parm_1.0.hpp"
//~ #include <string.h>
//~ #include <fstream>
//~ #include <iostream>
//~ #include <stdlib.h>
//~ #include <stdio.h>


#define z_axis_modify  2
#define deviation_y_center 10
#define deviation_x_center 10
#define start_nc 2
#define end_nc  13
#define diffusion_time 1000
#define name_parm7 "density_dis9a5.parm7" 
//~ #define name_nc "water_ion_graphene_10a5"

int main()
{
	typedef std::vector<double>::size_type index;
	std::cout << "program to calculate the diffusion constant"<< "\n"<< std::endl;
	
	//~ double water
	
    std::vector<std::vector<double>>  waterdistancesquare; 
    std::vector<double> aaa(3,0);
    std::vector<index>  total_water(diffusion_time,0);
    
    
    for(index i =0; i !=  diffusion_time; ++i)
    {
		waterdistancesquare.push_back(aaa);
		//~ std::cout << waterdistancesquare[i][2] <<std::endl;
	}			
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
	    std::cout<< "total water:  " << O_WAT_id.size() <<"\n" <<std::endl;
	    std::vector<index> O_WAT_IN_C_id;
	    	         
	   double C_z_coor_sum_1 = 0, C_z_coor_sum_2 = 0;
       double C_z_coor_average_1, C_z_coor_average_2;	
	   double C_x_coor_sum = 0, C_y_coor_sum = 0;
       double C_x_coor_center, C_y_coor_center;    
       
       double Z_UP,Z_DOWN,Y_UP,Y_DOWN,X_UP,X_DOWN;
     
       for(int frame =0; frame != total_frame; ++frame)
       //~ for(int frame =0; frame != 2       ;++frame)
       {
		   if((frame %1000) == 0)
		   {
		   std::cout << "frame_now :                                     " << frame <<std::endl;
	       }
		    if (frame==0)
		    {
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
	             //~ std::cout << "C_z_coor_average_1: " <<  C_z_coor_average_1 << std::endl;
	             //~ std::cout << "C_z_coor_average_2: " <<  C_z_coor_average_2 << std::endl;
	             //~ std::cout << "C_y_coor_center: " <<  C_y_coor_center << std::endl;
	             //~ std::cout << "C_x_coor_center: " <<  C_x_coor_center << std::endl;
	            
	             Z_UP = C_z_coor_average_2;
                 Z_DOWN = C_z_coor_average_1;
                 Y_UP = C_y_coor_center + deviation_y_center;
                 Y_DOWN = C_y_coor_center - deviation_y_center;
                 X_UP = C_x_coor_center + deviation_x_center ;
                 X_DOWN = C_x_coor_center - deviation_x_center ;
                 //~ std::cout << Z_UP <<std::endl;
                 //~ std::cout << Z_DOWN <<std::endl;      
                 //~ std::cout << Y_UP <<std::endl;
                 //~ std::cout << Y_DOWN <<std::endl;         
                 //~ std::cout << X_UP <<std::endl;  
                 //~ std::cout << X_DOWN<<std::endl;
	          }
	          
	         //~ std::vector<index> O_WAT_IN_C_id;
	         std::vector<double> O_coor,O_coor_initial;
	
             if( (frame %  diffusion_time)  == 0)
             //~ if(frame == 0)
             {
				 //~ std::cout << "frame % diffusion: "<< frame % diffusion_time << std::endl;
				 O_WAT_IN_C_id.clear();
				 for(index i = 0; i != O_WAT_id.size(); ++i)
				 {
					 O_coor = nc_data.atom_coordinate(frame, O_WAT_id[i]);
					 if( O_coor[2] < Z_UP && O_coor[2] > Z_DOWN && O_coor[0] < X_UP && O_coor[0] > X_DOWN && O_coor[1] < Y_UP && O_coor[1] > Y_DOWN)
					 {
						 O_WAT_IN_C_id.push_back(O_WAT_id[i]);
	                     //~ std::cout << O_WAT_id[i] << std::endl;
					  }
		          }
		          //~ std::cout<< " num_wat_in_c  "  << O_WAT_IN_C_id.size() << std::endl;
			  }
			  for(index i =0; i  != O_WAT_IN_C_id.size(); ++i)
			  //~ for(int i=0; i != 2; ++i)
			  {
				  //~ std::cout<< "water_id : " << O_WAT_IN_C_id[i] << std::endl;
				  //~ std::cout << "total_distance: "<< waterdistancesquare[frame % diffusion_time][0] << std::endl;
				  //~ std::cout << "coordinates_abs: " << nc_data.atom_coordinate(frame,O_WAT_IN_C_id[i])[0]<<std::endl;
				  //~ std::cout << "coordinates_ralative: " << nc_data.atom_coordinate(frame - frame % diffusion_time,O_WAT_IN_C_id[i])[0]<<std::endl;	
				  //~ std::cout << "coordinates: " << nc_data.atom_coordinate(frame,O_WAT_IN_C_id[i])[0]- nc_data.atom_coordinate(frame - frame % diffusion_time,O_WAT_IN_C_id[i])[0]<<std::endl;				  
				  //~ std::cout << "coordinates_square: "<< pow((nc_data.atom_coordinate(frame,O_WAT_IN_C_id[i])[0]- nc_data.atom_coordinate(frame - frame % diffusion_time,O_WAT_IN_C_id[i])[0]),2) <<std::endl;
				   waterdistancesquare[frame % diffusion_time][0] +=  pow((nc_data.atom_coordinate(frame,O_WAT_IN_C_id[i])[0]- nc_data.atom_coordinate(frame - frame % diffusion_time,O_WAT_IN_C_id[i])[0]),2);
				   waterdistancesquare[frame % diffusion_time][1] +=  pow((nc_data.atom_coordinate(frame,O_WAT_IN_C_id[i])[1]- nc_data.atom_coordinate(frame - frame % diffusion_time,O_WAT_IN_C_id[i])[1]),2);
				   waterdistancesquare[frame % diffusion_time][2] +=  pow((nc_data.atom_coordinate(frame,O_WAT_IN_C_id[i])[2]- nc_data.atom_coordinate(frame - frame % diffusion_time,O_WAT_IN_C_id[i])[2]),2);		
				   total_water[frame % diffusion_time] += 1;   		   
		      }
		      //~ std::cout << "waterdistancesquar:     "<<waterdistancesquare[frame % diffusion_time][0] << std::endl;	
	     }       // frame end     
	}      // nc end
	//~ for(int i = 0; i != diffusion_time; ++i)
	//~ {
		//~ std::cout << total_water[i] <<std::endl;
	//~ }
	std::ofstream outfile;
	outfile.open("data_diffusion");
	for(int i = 1; i != diffusion_time; ++i)
	{
		//~ std::cout << i * 0.25 << "    " << waterdistancesquare[i][0]/total_water[i]<< "    " << waterdistancesquare[i][1]/total_water[i]<<"    " << waterdistancesquare[i][0]/total_water[i]<<std::endl;
		//~ outfile << i * 0.25 << "    " << waterdistancesquare[i][0]/total_water[i]<< "    " << waterdistancesquare[i][1]/total_water[i]<<"    " << waterdistancesquare[i][0]/total_water[i]<<std::endl;		
		std::cout <<std::setw(10) << i * 0.25 << std::setw(10) << waterdistancesquare[i][0]/total_water[i]<< std::setw(10) << waterdistancesquare[i][1]/total_water[i]<<std::setw(10) << waterdistancesquare[i][2]/total_water[i]<<std::endl;
		outfile <<std::setw(10) << i * 0.25 << std::setw(10) << waterdistancesquare[i][0]/total_water[i]<< std::setw(10) << waterdistancesquare[i][1]/total_water[i]<<std::setw(10) << waterdistancesquare[i][2]/total_water[i]<<std::endl;
     }
     outfile.close();

	return 0;
}
