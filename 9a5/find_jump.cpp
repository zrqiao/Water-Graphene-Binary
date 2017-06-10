#include "amber_netcdf.hpp"
#include "amber_parm_1.0.hpp"
//~ #include <string.h>
//~ #include <fstream>
//~ #include <iostream>
//~ #include <stdlib.h>
//~ #include <stdio.h>



#define deviation_y_center 10
#define deviation_x_center 10
#define start_nc 2
#define end_nc  29
#define name_parm7 "density_dis9a5.parm7"
#define jump_time 1000
#define dt 4
//~ #define name_nc "water_ion_graphene_10a5"
typedef std::vector<double>::size_type index;


int main() {
    std::vector<int*> jump_coor;
    std::ofstream outfile;
    outfile.open("find_jump_coor");
	std::cout << "program to calculate the molecules that have jumped" << "\n" << std::endl;


	long total_wat = 0;

	for (int nc = start_nc; nc != end_nc + 1; ++nc) {
		char name_nc[64];
		sprintf(name_nc, "density_dis9a5_%d.nc", nc);
		amber_parm parm_name(name_parm7);
		nctraj nc_data(name_nc);
		printf("nc_file = %d\n", nc);
		int total_frame = nc_data.frames_number();
		std::cout << "total frame: " << total_frame << std::endl;
		double C_z_coor_sum_1 = 0, C_z_coor_sum_2 = 0;
		double C_z_coor_average_1, C_z_coor_average_2;
		double C_x_coor_sum = 0, C_y_coor_sum = 0;
		double C_x_coor_center, C_y_coor_center;
		double Z_UP, Z_DOWN, Y_UP, Y_DOWN, X_UP, X_DOWN, MID_DOWN, MID_UP;


		std::vector<index> O_WAT_id = parm_name.id_by_type("OW");
		std::vector<index> O_WAT_IN_C_id;
		std::vector<index> WAT_JUMP;
		std::vector<double> O_coor, O_coor_initial;
        for (index C_index = 0; C_index != 1400; ++C_index) {
            C_z_coor_sum_1 += nc_data.atom_coordinate(0, C_index)[2];
            C_z_coor_sum_2 += nc_data.atom_coordinate(0, (1400 + C_index))[2];
            C_y_coor_sum += nc_data.atom_coordinate(0, C_index)[1];
            C_x_coor_sum += nc_data.atom_coordinate(0, C_index)[0];
        }
        C_z_coor_average_1 = C_z_coor_sum_1 / 1400;
        C_z_coor_average_2 = C_z_coor_sum_2 / 1400;
        C_y_coor_center = C_y_coor_sum / 1400;
        C_x_coor_center = C_x_coor_sum / 1400;

        Z_UP = C_z_coor_average_2;
        Z_DOWN = C_z_coor_average_1;
        Y_UP = C_y_coor_center + deviation_y_center;
        Y_DOWN = C_y_coor_center - deviation_y_center;
        X_UP = C_x_coor_center + deviation_x_center;
        X_DOWN = C_x_coor_center - deviation_x_center;
        MID_DOWN = Z_DOWN + (Z_UP - Z_DOWN)*7/15;
        MID_UP = Z_UP - (Z_UP - Z_DOWN)*7/15;
        for (index frame0=0;frame0<total_frame-jump_time+1;frame0+=(jump_time)) {
            int jump_count=0;
            //std::cout<<"frame now: "<<frame0<<std::endl;

            O_WAT_IN_C_id.clear();
            for (index i = 0; i != O_WAT_id.size(); ++i) {
                O_coor = nc_data.atom_coordinate(0, O_WAT_id[i]);
                if (O_coor[2] < Z_UP && O_coor[2] > Z_DOWN && O_coor[0] < X_UP && O_coor[0] > X_DOWN &&
                    O_coor[1] < Y_UP && O_coor[1] > Y_DOWN) {
                    O_WAT_IN_C_id.push_back(O_WAT_id[i]);
                }
            }

            total_wat += O_WAT_IN_C_id.size();

            int layerflag_init=0;
            int layer=0;
            double z_coor;
            //std::cout<<"calculate the target atoms"<<"\n" << std::endl;

            for (index i = 0; i != O_WAT_IN_C_id.size(); i++) {
                //std::cout<<O_WAT_IN_C_id[i]<<std::endl;
                for (index frame = frame0; frame!= (frame0+jump_time); frame+=dt) {

                    z_coor = nc_data.atom_coordinate(frame, O_WAT_IN_C_id[i])[2];

                    if (Z_DOWN <= z_coor&&z_coor < MID_DOWN) {
                        layer= -1;

                    } else if (MID_DOWN <=z_coor&& z_coor < MID_UP) {
                        layer= 0;

                    } else if (MID_UP <= z_coor&& z_coor <= Z_UP) {
                        layer= 1;

                    }
                    else{}
                    //std::cout<<layer<<' '<<MID_DOWN<<std::endl;
                    if (layer!= 0) { //初始位置不能在中心
                        layerflag_init = layer;

                        break;
                    }
                }
                for (index frame = frame0; frame != frame0+jump_time; frame+=dt) {

                    z_coor = nc_data.atom_coordinate(frame, O_WAT_IN_C_id[i])[2];
                    if (Z_DOWN <= z_coor&& z_coor< MID_DOWN) {
                        layer= -1;
                    } else if (MID_DOWN <= z_coor&& z_coor< MID_UP) {
                        layer= 0;
                    } else if (MID_UP <= z_coor&& z_coor<= Z_UP) {
                        layer= 1;
                    }

                    if (layer  == -1 && layerflag_init == 1) {
                        O_coor = nc_data.atom_coordinate(frame, O_WAT_IN_C_id[i]);
                        if (O_coor[2] < Z_UP && O_coor[2] > Z_DOWN && O_coor[0] < X_UP && O_coor[0] > X_DOWN &&
                            O_coor[1] < Y_UP && O_coor[1] > Y_DOWN) {
                            jump_count++;
                            outfile<<nc*total_frame+frame<<std::setw(10)<<O_coor[0]<<std::setw(10)<<O_coor[1]<<std::endl;
                            std::cout<<nc*total_frame+frame<<std::setw(10)<<O_coor[0]<<std::setw(10)<<O_coor[1]<<std::endl;

                        }

                    }
                    if (layer*layerflag_init==-1){
                        layerflag_init = layer;
                    }
                    else{
                        //std::cout << frame<<' '<<0<< "\n" << std::endl;
                        //outfile << frame<<' '<<0<< "\n" << std::endl;
                    }

                }

            }
            std::cout << frame0<<' '<<jump_count << std::endl;
            //outfile << frame0<<' '<<jump_count << std::endl;
        }

		       // frame start


	}      // nc end
    outfile.close();

	return 0;
}
