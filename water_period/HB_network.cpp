//
// Created by utena on 17-10-15.
//
#include "amber_netcdf.hpp"
#include "amber_parm_1.0.hpp"
#include "vector_calc.h"
#define start_nc 37
#define end_nc  100
#define select_boundary 3.55
#define hbond_cutoff_up 3.5
#define hbond_cutoff_down 2.0
#define hbond_cutoff_angle_up 180.0
#define hbond_cutoff_angle_down 135.0
#define period_boundary
#define x_points 15
#define y_points 15
#define dt 4
#define dr 0.01
#define name_parm7 "density_dis6a5.parm7"
typedef std::vector<double>::size_type index;
int main() {
    std::vector<double> OH1_vec(3,0),OH2_vec(3,0),Dip_vec(3,0),Nor_vec(3,0),O1H1_vec(3,0),O2H2_vec(3,0);
    char name_nc[64];
    sprintf(name_nc, "nc/density_dis6a5_%d.nc", start_nc);
    amber_parm parm_name(name_parm7);
    nctraj nc_data(name_nc);
    std::vector<index> O_WAT_id = parm_name.id_by_type("OW");
    std::vector<index> O_WAT_IN_C_id, select_wat_id;
    double C_x_coor_sum = 0, C_y_coor_sum = 0;
    double C_x_coor_center, C_y_coor_center;
    std::vector<double> O_coor, O1_coor,O2_coor, H1O_coor, H2O_coor,H1_coor, H2_coor;
    double angleH1,angleH2,angleO1H1,angleO2H2,angleDip,angleNor;
    for (int nc=start_nc;nc<=end_nc;nc++) {
        amber_parm parm_nam(name_parm7);
        char name_nc[64];
        sprintf(name_nc, "nc/density_dis6a5_%d.nc",nc);
        nctraj data_nc(name_nc);
        for (index C_index = 0; C_index != 1400; ++C_index) {
            C_y_coor_sum += data_nc.atom_coordinate(0, C_index)[1];
            C_x_coor_sum += data_nc.atom_coordinate(0, C_index)[0];
        }
    }
    C_z_coor_average_1 = C_z_coor_sum_1/(1400*(end_nc-start_nc+1));
    C_z_coor_average_2 = C_z_coor_sum_2/(1400*(end_nc-start_nc+1));
    C_y_coor_center = C_y_coor_sum/(1400*(end_nc-start_nc+1));
    C_x_coor_center = C_x_coor_sum/(1400*(end_nc-start_nc+1));

    Y_UP = C_y_coor_center + deviation_y_center;
    Y_DOWN = C_y_coor_center - deviation_y_center;
    X_UP = C_x_coor_center + deviation_x_center ;
    X_DOWN = C_x_coor_center - deviation_x_center ;

    std::cout << "X_UP: " << X_UP <<std::endl;
    std::cout << "X_DOWN: " << X_DOWN <<std::endl;
    std::cout << "Y_UP: " << Y_UP <<std::endl;
    std::cout << "Y_DOWN: " << Y_DOWN <<std::endl;
    int Target_ID,nc,frame,HB_ID1,HB_ID2;
    for (index frame=0;frame<nc_data.frames_number();frame+=dt) {
        std::vector<std::vector<std::vector<double>>> xypoints(x_points+1,std::vector<std::vector<double>>(y_points+1,new std::vector<double>));
        if (frame%100==0) std::cout<<frame<<std::endl;
        O_WAT_IN_C_id.clear();select_wat_id.clear();
        for(index i = 0; i != O_WAT_id.size(); ++i) {//打格点
            O_coor = nc_data.atom_coordinate(frame, O_WAT_id[i]);
            int x_location=floor(O_coor[0]/select_boundary);
            int y_location=floor(O_coor[1]/select_boundary);
            xypoints[x_location][y_location].push_back(O_WAT_id[i]);
            xypoints[x_location+1][y_location].push_back(O_WAT_id[i]);
            xypoints[x_location][y_location+1].push_back(O_WAT_id[i]);
            xypoints[x_location+1][y_location+1].push_back(O_WAT_id[i]);
        }
        for (int i=0;i<x_points;i++) {
            xypoints[i][0].insert(xypoints[i][0].end(),xypoints[i][y_points].begin(),xypoints[i][ypoints].end());
        }
        for (int j=0;j<y_points;j++) {
            xypoints[0][j].insert(xypoints[0][j].end(),xypoints[x_points][j].begin(),xypoints[x_points][j].end());
        }
        for (int i=0;i<x_points;i++) {
            for (int j=0;j<y_points;j++) {
                for (vector<double>::iterator iter=xypoints[i][j].begin();iter!=xypoints[i][j].end();iter++){
                    int wat1=*iter;
                    for (vector<double>::iterator iter2=xypoints[i][j].begin();iter2!=xypoints[i][j].end();iter2++){
                        int wat2=*iter2;
                        O1_coor = nc_data.atom_coordinate(frame, wat1);
                        O2_coor = nc_data.atom_coordinate(frame, wat2);
                        H1_coor = nc_data.atom_coordinate(frame, wat1 + 1);
                        H2_coor = nc_data.atom_coordinate(frame, wat1 + 2);
                        if (judge_H_Bond(O1_coor, H1_coor, O2_coor)){
                            //输出
                        }
                        if (judge_H_Bond(O1_coor, H2_coor, O2_coor)){

                        }
                    }
                }
            }
        }
    }
}

