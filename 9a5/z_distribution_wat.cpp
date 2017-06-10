#include "amber_netcdf.hpp"
#include "amber_parm_1.0.hpp"
#include "vector_calc.h"


#define z_axis_modify  2
#define deviation_y_center 20
#define deviation_x_center 20
#define start_nc 0
#define end_nc  0
#define z_axis_num 100
#define name_parm7 "density_dis16a5.parm7"
//~ #define name_nc "water_ion_graphene_10a5"

int main()
{
	typedef std::vector<double>::size_type index;
	std::cout << "program to calculate the Z_axis distribution"<< "\n"<< std::endl;
	std::cout << "calculate XYZ_limit" << std::endl;
	
	std::vector<int> z_axis_dis(z_axis_num,0);
	
    double C_z_coor_sum_1 = 0, C_z_coor_sum_2 = 0;
    double C_z_coor_average_1, C_z_coor_average_2;	
	double C_x_coor_sum = 0, C_y_coor_sum = 0;
    double C_x_coor_center, C_y_coor_center;     
    double Z_UP,Z_DOWN,Y_UP,Y_DOWN,X_UP,X_DOWN;
	
	amber_parm parm_nam(name_parm7);
	nctraj data_nc("density_dis16a5_0.nc");

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
    
	for(int nc = start_nc; nc  != end_nc+1; ++nc)
	{
        char name_nc[64];
        sprintf(name_nc, "density_dis16a5_%d.nc",nc);
		amber_parm parm_name(name_parm7);
        nctraj nc_data(name_nc);
        printf("nc_file = %d\n",nc);
        int total_frame = nc_data.frames_number();
        std::vector<index> O_WAT_id = parm_name.id_by_type("OW");	
        int frame = 0;
       std::vector<double> o_coor;
       double step;
       while(frame 	<     total_frame)
       {
		    for(index i = 0; i != O_WAT_id.size(); ++i)
		    {
		        o_coor = nc_data.atom_coordinate(frame, O_WAT_id[i]);
				if( o_coor[2] < Z_UP && o_coor[2] > Z_DOWN && o_coor[0] < X_UP && o_coor[0] > X_DOWN && o_coor[1] < Y_UP && o_coor[1] > Y_DOWN)
				{
					step = (Z_UP-Z_DOWN)/z_axis_num;
				    z_axis_dis[floor((o_coor[2] - Z_DOWN)/step)] += 1;
				}
			}
		   frame += 10;
	   }
    }
	
	

	
	std::ofstream outfile;

	
	outfile.open("z_axis_distribution");

	for(index i = 0; i != z_axis_dis.size(); ++i)
	{
		outfile <<std::setw(15) <<i<<std::setw(15)<<z_axis_dis[i] << std::endl;
		std::cout <<std::setw(15) <<i<<std::setw(15)<<z_axis_dis[i] << std::endl;
     }
     outfile.close();      
	return 0;
}






