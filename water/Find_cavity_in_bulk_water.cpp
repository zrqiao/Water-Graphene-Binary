//
// Created by utena on 17-6-29.
//

#include "amber_netcdf.hpp"
#include "amber_parm_1.0.hpp"
#include<cmath>
#include<cstdlib>
#include<cstring>
#include<cfloat>
#define z_axis_modify  2
#define y_length 20
#define x_length 20
#define start_nc 1
#define end_nc  1
#define dt 1000
#define name_parm7 "water_tip3p.prmtop"
#define psi 2.4
#define cavity_cut 0.008
#define dx 0.4
#define cut_off 24
#define max_size 200
#define max_life 10000
#define rapid_change_constant 5
static int x_points = 50;
static int y_points = 50;

double density(double x,double y,double z, double x_W, double y_W,double z_W){
    double r= sqrt(pow((x-x_W),2)+pow(y-y_W,2)+pow((z-z_W),2));
    return (1/(2*M_PI*pow(psi,2)))*exp(-2*(pow((r/psi),2)));
}

int calc_cavity_size(int cavity_num,int x,int y,int z,int &last_step_size, int*** &judge_matrix, int*** &calc_matrix,std::vector<int*> &temp_cavity_coordinate,int x_points, int y_points){
    int *coordinate_vector=new int [3];
    if (calc_matrix[x][y][z]==0) {
        last_step_size += 1;
        calc_matrix[x][y][z] = 1;
        for (int i=0;i<3;i++){
            coordinate_vector[0] = x;
            coordinate_vector[1] = y;
            coordinate_vector[2] = z;
            coordinate_vector[i]--;
            /*if (coordinate_vector[i]<0){
                coordinate_vector+=x_points;
            }*/
            if (0<=coordinate_vector[i]&&coordinate_vector[i]<x_points) {
                if (judge_matrix[coordinate_vector[0]][coordinate_vector[1]][coordinate_vector[2]] == 1 &&
                    calc_matrix[coordinate_vector[0]][coordinate_vector[1]][coordinate_vector[2]] == 0) {
                    temp_cavity_coordinate.push_back(coordinate_vector);
                    calc_cavity_size(cavity_num, coordinate_vector[0], coordinate_vector[1], coordinate_vector[2],
                                     last_step_size, judge_matrix, calc_matrix, temp_cavity_coordinate,
                                     x_points, y_points);
                }
            }
            coordinate_vector[0] = x;
            coordinate_vector[1] = y;
            coordinate_vector[2] = z;
            coordinate_vector[i]++;
            /*if (coordinate_vector[i]>=x_points){
                coordinate_vector-=x_points;
            }*/
            if (0<=coordinate_vector[i]&&coordinate_vector[i]<x_points) {
                if (judge_matrix[coordinate_vector[0]][coordinate_vector[1]][coordinate_vector[2]] == 1 &&
                    calc_matrix[coordinate_vector[0]][coordinate_vector[1]][coordinate_vector[2]] == 0) {
                    temp_cavity_coordinate.push_back(coordinate_vector);
                    calc_cavity_size(cavity_num, coordinate_vector[0], coordinate_vector[1], coordinate_vector[2],
                                     last_step_size, judge_matrix, calc_matrix, temp_cavity_coordinate,
                                     x_points, y_points);

                }
            }
        }
    }
    return last_step_size;
}
int main() {

    std::ofstream outfile2;
    typedef std::vector<double>::size_type index;
    std::cout << "program to find the cavities" << "\n" << std::endl;

    double max_cavity_index=0;
    std::vector<double> cavity_quantities_by_size(max_size,0);
    std::vector<double> cavity_quantities_by_life(max_life,0);

    std::vector<int> active_cavity_index;
    std::vector<int> new_active_index;
    std::vector<std::vector<int*>> last_step_cavity_coordinate_by_index;
    static std::vector<std::vector<int>> cavity_size_by_time;
    static std::vector<int> cavity_frame;
    static int cavity_num=0;
    static double ***density_matrix = new double **[x_points];
    static int ***judge_matrix = new int **[x_points];
    static int ***calc_matrix = new int **[x_points];
    for (int i = 0; i < x_points; i++) {
        density_matrix[i] = new double*[x_points];
        judge_matrix[i]=new int *[x_points];
        calc_matrix[i] = new int* [y_points];
        for (int j=0; j<x_points;j++){
            density_matrix[i][j]=new double [x_points];
            judge_matrix[i][j]=new int [x_points];
            calc_matrix[i][j] = new int [y_points];
        }
    }

    for (int i = 0; i < x_points; i++) {

        for (int j=0; j<x_points;j++){
            density_matrix[i][j]=new double [x_points];
        }
    }

    for (int i = 0; i < x_points; i++) {

    }
    static std::vector<double> cavity_by_index;
    std::vector<int*> temp_coordinates;
    for (int nc = start_nc; nc != end_nc + 1; ++nc) {
        char name_nc[64];
        char name_out2[64];
        sprintf(name_nc, "water_tip3p_%d.nc", nc);
        sprintf(name_out2,"cavity_size_by_time_%d",nc);
        amber_parm parm_name(name_parm7);
        nctraj nc_data(name_nc);
        printf("nc_file = %d\n", nc);
        int total_frame = nc_data.frames_number();
        std::cout << "total frame: " << total_frame << std::endl;
        outfile2.open(name_out2);
        std::vector<index> O_WAT_id = parm_name.id_by_type("OW");
        std::cout << "total water:  " << O_WAT_id.size() << "\n" << std::endl;
        std::vector<index> O_WAT_IN_C_id_upperlayer;
        std::vector<index> O_WAT_IN_C_id_lowerlayer;

        double Z_UP=30, Z_DOWN=10, Y_UP=30, Y_DOWN=10, X_UP=30, X_DOWN=10;

        for (int frame = 0; frame != total_frame; frame+=dt)
            //~ for(int frame =0; frame != 2       ;++frame)
        {
            if ((frame % 1000) == 0) {
                std::cout << "frame_now : " << frame << std::endl;
            }
            std::vector<double> O_coor,O_coor_initial;
            if(frame %dt==0 && frame!=total_frame) {
                O_WAT_IN_C_id_lowerlayer.clear();
                O_WAT_IN_C_id_upperlayer.clear();
                for (int i = 0; i < x_points; i++) {
                    for (int j = 0; j < x_points; j++) {
                        for(int k=0; k < x_points; k++) {
                            density_matrix[i][j][k] = 0;
                            judge_matrix[i][j][k] = 0;
                        }
                    }
                }
                double temp_O_X, temp_O_Y, temp_O_Z;
                for(index m =0; m  != O_WAT_id.size(); ++m) {
                    temp_O_X=0;
                    temp_O_Y=0;
                    temp_O_Z=0;
                    for (int f=frame;f<frame+dt;f++){
                        temp_O_X+=(nc_data.atom_coordinate(f,O_WAT_id[m])[0]/dt);
                        temp_O_Y+=(nc_data.atom_coordinate(f,O_WAT_id[m])[1]/dt);
                        temp_O_Z+=(nc_data.atom_coordinate(f,O_WAT_id[m])[2]/dt);
                    }
                    for (int i=0; i<x_points;i++){
                        for(int j=0;j<x_points;j++){
                            for(int k=0;k<x_points;k++){
                                density_matrix[i][j][k]+=density(X_DOWN+i*dx,Y_DOWN+j*dx,Z_DOWN+k*dx,temp_O_X,temp_O_Y,temp_O_Z);//计算密度分布函数
                            }
                        }

                    }
                    //std::cout<<m<<std::endl;
                }
                std::cout<<"step: "<<frame/dt<<'\n'<<std::endl;
                for (int i=0; i<x_points;i++){
                    for(int j=0;j<x_points;j++){
                        for(int k=0;k<x_points;k++) {
                            if (density_matrix[i][j][k] < cavity_cut) {
                                judge_matrix[i][j][k] = 1;
                            } else {
                                judge_matrix[i][j][k] = 0;
                            }
                        }
                        //通过密度矩阵判断是否为空穴
                        //std::cout<<judge_matrix[i][j][0]<<' ';
                    }
                    //std::cout<<std::endl;
                }
                temp_coordinates.clear();
                bool is_active= false;
                new_active_index.clear();
                int cavity_size=0;
                for (int i=0; i<x_points;i++){
                    for(int j=0;j<x_points;j++){
                        for (int k=0;k<x_points;k++){
                            if (judge_matrix[i][j][k]==1&& calc_matrix[i][j][k]==0) {
                                cavity_size = 0;
                                temp_coordinates.clear();
                                calc_cavity_size(cavity_num, i, j, k, cavity_size, judge_matrix, calc_matrix, temp_coordinates, x_points, y_points);//
                                /*if (cavity_size>500){
                                    std::cout<<frame<<std::setw(10)<<cavity_size<<std::endl;
                                    for (int m=0;m<O_WAT_id.size();m++){
                                        for (int c=0;c<temp_coordinates.size();c++){
                                            if (abs(nc_data.atom_coordinate(frame,O_WAT_id[m])[0]-(X_DOWN+dx*temp_coordinates[c][0]))<=dx &&
                                                    abs(nc_data.atom_coordinate(frame,O_WAT_id[m])[1]-(Y_DOWN+dx*temp_coordinates[c][1]))<=dx &&
                                                        abs(nc_data.atom_coordinate(frame,O_WAT_id[m])[2]-(Z_DOWN+dx*temp_coordinates[c][2]))<=dx){
                                                std::cout<<O_WAT_id[m]<<std::setw(5);
                                            }
                                        }
                                    }
                                    std:: cout<<std::endl;
                                }*/
                                cavity_quantities_by_size[cavity_size]++;
                                cavity_num++;
                            }
                            //std::cout<<cavity_size<<std::endl;
                        }
                    }
                }
                std::ofstream outfile;
                char out_name[64];
                sprintf(out_name,"cavity_size_distribution_bulk_%dfs.dat",dt*4);
                outfile.open(out_name);
                for (index i=0; i<max_size;i++){
                    std::cout<<i*pow(dx,3)<<std::setw(15)<<cavity_quantities_by_size[i]*dt/(frame+dt)<<std::endl;
                    outfile<<i*pow(dx,3)<<std::setw(15)<<cavity_quantities_by_size[i]*dt/(frame+dt)<<std::endl;
                }
                outfile.close();
                //std::cout<<last_step_cavity_coordinate_by_index.size()<<std::endl;
                /*for (int i=0;i<active_cavity_index.size();i++){//计算空穴覆盖情况
                    is_active=false;
                    int temp_x=0,temp_y=0,temp_z=0;
                    for (;!last_step_cavity_coordinate_by_index[active_cavity_index[i]].empty();){
                        int *temp_xy = last_step_cavity_coordinate_by_index[active_cavity_index[i]].back();
                        last_step_cavity_coordinate_by_index[active_cavity_index[i]].pop_back();
                        if (judge_matrix[temp_xy[0]][temp_xy[1]][temp_xy[2]] == 1) {
                            is_active = true;//版本1 有重合就是同一个
                            temp_x = temp_xy[0];
                            temp_y = temp_xy[1];
                            temp_z = temp_xy[2];
                            //std::cout<<temp_x<<' '<<temp_y<<' ';
                        }

                    }
                    last_step_cavity_coordinate_by_index[active_cavity_index[i]].clear();
                    cavity_size=0;
                    if (is_active){

                        calc_cavity_size(active_cavity_index[i],temp_x,temp_y,temp_z,cavity_size,judge_matrix,calc_matrix,last_step_cavity_coordinate_by_index[active_cavity_index[i]],x_points,x_points);

                    }
                    new_active_index.push_back(active_cavity_index[i]);
                    if (abs(cavity_size-cavity_size_by_time[active_cavity_index[i]].back())*pow(dx,2)>rapid_change_constant){
                        outfile<<total_frame*nc+frame<<std::setw(7)<<active_cavity_index[i]<<std::setw(7)<<pow(dx,2)*abs(cavity_size-cavity_size_by_time[active_cavity_index[i]].back())<<std::setw(7)<<temp_x<<std::setw(7)<<temp_y<<std::endl;
                    }
                    cavity_size_by_time[active_cavity_index[i]].push_back(cavity_size);
                    if (cavity_size==0){
                        cavity_quantities_by_life[cavity_size_by_time[active_cavity_index[i]].size()]++;
                        //for (index j=0;j<cavity_size_by_time[active_cavity_index[i]].size();j++){
                        //  outfile2<<cavity_size_by_time[active_cavity_index[i]][j]*pow(dx,2)<<std::setw(10);
                        //}
                        cavity_size_by_time[active_cavity_index[i]].clear();
                        outfile2<<frame<<std::endl;
                    }
                }
                active_cavity_index.clear();
                for (int m=0;m<new_active_index.size();m++){//更新活空穴
                    active_cavity_index.push_back(new_active_index[m]);
                }
                for (int i=0; i<x_points;i++){//深搜统计新的空穴大小和数目
                    for(int j=0;j<x_points;j++){
                        for (int k=0;k<x_points;k++){
                        if (judge_matrix[i][j][k]==1&& calc_matrix[i][j]==0) {
                            cavity_size = 0;
                            temp_coordinates.clear();
                            calc_cavity_size(cavity_num, i, j, k, cavity_size, judge_matrix, calc_matrix, temp_coordinates, x_points, y_points);//
                            //temp_coordinates.clear();//完全忽略面积信息，内存占用极大时使用
                            //std::cout<<temp_coordinates.size()<<std::endl;
                            last_step_cavity_coordinate_by_index.push_back(temp_coordinates);
                            std::vector<int> init_cavity_size_vector;
                            init_cavity_size_vector.push_back(cavity_size);
                            cavity_size_by_time.push_back(init_cavity_size_vector);
                            cavity_frame.push_back(frame);
                            active_cavity_index.push_back(cavity_num);
                            cavity_quantities_by_size[cavity_size]++;
                            cavity_num++;
                        }
                            //std::cout<<cavity_size<<std::endl;
                        }
                    }
                }*/
                //std::cout<<frame <<' '<<large_cavity_num<<std::endl;

            }
        }
        outfile2.close();
    }
    std::ofstream outfile1;


    /*for (int i=0; i<x_points;i++) {
        for (int j = 0; j < y_points; j++) {
            std::cout << i * dx << std::setw(10) << j * dx << std::setw(10) << sum_matrix[i][j] << std::endl;
        }
    }*/
    /*outfile1.open("cavity_lifetime");

    for (index i=0;i<cavity_size_by_time.size();i++){


    }
    for (index i=0; i<max_life;i++) {
        outfile1 << i << std::setw(10) << cavity_quantities_by_life[i] << std::endl;
    }
    outfile1.close();*/


    return 0;
}
