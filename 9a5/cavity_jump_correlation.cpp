//
// Created by utena on 17-5-28.
//

//
// Created by utena on 17-5-8.
//
#include "enthaply/amber_netcdf.hpp"
#include "enthaply/amber_parm_1.0.hpp"
#include<cmath>
#include<cstdlib>
#include<cstring>
#include<cfloat>
#include <algorithm>

#define z_axis_modify  2
#define deviation_y_center 10
#define deviation_x_center 10
#define start_nc 2
#define end_nc  3
#define dt 1
#define name_parm7 "density_dis9a5.parm7"
#define psi 2.4
#define cavity_cut 0.016
#define dx 0.2
#define cut_off 24
#define max_size 1000
#define max_life 1000

double density(double x,double y, double x_W, double y_W){
    double r= sqrt(pow((x-x_W),2)+pow((y-y_W),2));
    return (1/(2*M_PI*pow(psi,2)))*exp(-2*(pow((r/psi),2)));
}
int calc_cavity_size(int cavity_num,int x,int y,int &last_step_size, int **judge_matrix, int **calc_matrix,std::vector<int*> &temp_cavity_coordinate,int x_points, int y_points){
    int *coordinate_vector=new int [2];
    if (calc_matrix[x][y]==0){
        last_step_size+=1;
        calc_matrix[x][y]=cavity_num+1;
        if (x-1>=0){
            if (judge_matrix[x-1][y]==1 && calc_matrix[x-1][y]==0){
                coordinate_vector[0]=x;
                coordinate_vector[1]=y;
                temp_cavity_coordinate.push_back(coordinate_vector);
                calc_cavity_size(cavity_num,x-1,y,last_step_size,judge_matrix,calc_matrix,temp_cavity_coordinate,x_points,y_points);
            }
        }
        if (y-1>=0){
            if (judge_matrix[x][y-1]==1 && calc_matrix[x][y-1]==0){
                coordinate_vector[0]=x;
                coordinate_vector[1]=y;
                temp_cavity_coordinate.push_back(coordinate_vector);
                calc_cavity_size(cavity_num,x,y-1,last_step_size,judge_matrix,calc_matrix,temp_cavity_coordinate,x_points,y_points);
            }
        }
        if (x+1<x_points){
            if (judge_matrix[x+1][y]==1 && calc_matrix[x+1][y]==0){
                coordinate_vector[0]=x;
                coordinate_vector[1]=y;
                temp_cavity_coordinate.push_back(coordinate_vector);
                calc_cavity_size(cavity_num,x+1,y,last_step_size,judge_matrix,calc_matrix,temp_cavity_coordinate,x_points,y_points);
            }
        }
        if (y+1<y_points){
            if (judge_matrix[x][y+1]==1 && calc_matrix[x][y+1]==0){
                coordinate_vector[0]=x;
                coordinate_vector[1]=y;
                temp_cavity_coordinate.push_back(coordinate_vector);
                calc_cavity_size(cavity_num,x,y+1,last_step_size,judge_matrix,calc_matrix,temp_cavity_coordinate,x_points,y_points);
            }
        }
    }
    return last_step_size;
}
int main() {
    int target_jump=0;
    int nt_jump=0;
    std::ifstream infile;
    //std::ifstream infile1;
    std::ofstream outfile;
    std::ofstream outfile2;
    infile.open("find_jump_coor");
    outfile.open("cavity_size_change");
    outfile2.open("size_change_distribution");
    //infile1.open("Rapid_size_change");
    typedef std::vector<double>::size_type index;
    static std::vector<double> cavity_frames;
    std::vector<double> density_distribution(500,0);
    std::vector<int> size_change_distribution(3000,0);
    std::cout << "program to derive correlations between jump& cavity" << "\n" << std::endl;
    int x_points = deviation_x_center * 2 / dx+1;
    int y_points = deviation_y_center * 2 / dx+1;
    static std::vector<int> cavity_frame;
    static double **density_matrix = new double *[x_points];
    static std::vector<int*> temp_coordinates;
    for (int i = 0; i < x_points; i++) {
        density_matrix[i] = new double[y_points];
    }
    static int **judge_matrix = new int *[x_points];
    for (int i = 0; i < x_points; i++) {
        judge_matrix[i] = new int[y_points];
    }
    static int **calc_matrix = new int *[x_points];
    for (int i = 0; i < x_points; i++) {
        calc_matrix[i] = new int[y_points];
    }
    double C_z_coor_sum_1 = 0, C_z_coor_sum_2 = 0;
    double C_z_coor_average_1, C_z_coor_average_2;
    double C_x_coor_sum = 0, C_y_coor_sum = 0;
    double C_x_coor_center, C_y_coor_center;

    double Z_UP, Z_DOWN, Y_UP, Y_DOWN, X_UP, X_DOWN;
    index frame=0;
    char name_nc[64];
    sprintf(name_nc, "density_dis9a5_%d.nc", 10);
    amber_parm parm_name(name_parm7);
    nctraj nc_data(name_nc);
    for (index C_index = 0; C_index != 1400; ++C_index) {
        C_z_coor_sum_1 += nc_data.atom_coordinate(frame, C_index)[2];
        C_z_coor_sum_2 += nc_data.atom_coordinate(frame, (1400 + C_index))[2];
        C_y_coor_sum += nc_data.atom_coordinate(frame, C_index)[1];
        C_x_coor_sum += nc_data.atom_coordinate(frame, C_index)[0];
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


    std::vector<double> O_coor,O_coor_initial;
    /*while (!infile1.fail()){
        double temp[5];
        infile1>>temp[0]>>temp[1]>>temp[2]>>temp[3]>>temp[4];
        cavity_frames.push_back(temp[0]);
    }*/
    while (!infile.fail()){
        double num[3];
        infile>>num[0]>>num[1]>>num[2];

        int frame_r=num[0];
        int nc=frame_r/10000;
        frame=frame_r-10000*nc-dt;
        //std::cout<<frame_r<<std::endl;
        int x_r=floor((num[1]-X_DOWN)/dx);
        int y_r=floor((num[2]-Y_DOWN)/dx);
        sprintf(name_nc, "density_dis9a5_%d.nc", nc);
        amber_parm parm_name(name_parm7);
        nctraj nc_data(name_nc);
        std::vector<index> O_WAT_id = parm_name.id_by_type("OW");
        std::vector<index> O_WAT_IN_C_id_upperlayer;
        std::vector<index> O_WAT_IN_C_id_lowerlayer;
        O_WAT_IN_C_id_lowerlayer.clear();
        O_WAT_IN_C_id_upperlayer.clear();
        int judge=0;
        /*for (int i=0;i<cavity_frames.size();i++){
            if (num[0]==cavity_frames[i]){
                judge=1;
                break;
            }
        }
        if (judge==1){
            target_jump++;
            std::cout<<frame_r<<' '<<x_r<<' '<<y_r<<std::endl;
        }
        else{
            nt_jump++;
        }*/
        for(index i = 0; i != O_WAT_id.size(); ++i)
        {
            O_coor = nc_data.atom_coordinate(frame, O_WAT_id[i]);
            if( O_coor[2] < Z_DOWN + (Z_UP - Z_DOWN)*7/15 && O_coor[2] >= Z_DOWN && O_coor[0] < X_UP && O_coor[0] > X_DOWN && O_coor[1] < Y_UP && O_coor[1] > Y_DOWN)//做修改，将水分为两层
            {
                O_WAT_IN_C_id_lowerlayer.push_back(O_WAT_id[i]);
                //~ std::cout << O_WAT_id[i] << std::endl;
            }
            if( O_coor[2] <Z_UP && O_coor[2] >=  Z_UP - (Z_UP - Z_DOWN)*7/15  && O_coor[0] < X_UP && O_coor[0] > X_DOWN && O_coor[1] < Y_UP && O_coor[1] > Y_DOWN)//做修改，将水分为两层
            {
                O_WAT_IN_C_id_upperlayer.push_back(O_WAT_id[i]);
                //~ std::cout << O_WAT_id[i] << std::endl;
            }
        }

        double temp_O_X; double temp_O_Y;
        for (int i = 0; i < x_points; i++) {
            for (int j = 0; j < y_points; j++) {
                density_matrix[i][j] = 0;
                judge_matrix[i][j] = 0;
                calc_matrix[i][j]=0;
                //初始化
            }
        }
        for(index m =0; m  != O_WAT_IN_C_id_lowerlayer.size(); ++m) {
            temp_O_X=0;
            temp_O_Y=0;
            /*for (int f=frame;f<frame+dt;f++){
                temp_O_X+=(nc_data.atom_coordinate(frame,O_WAT_IN_C_id_lowerlayer[m])[0]/dt);
                //std::cout<<temp_O_X<<'\n'<<std::endl;
                temp_O_Y+=(nc_data.atom_coordinate(frame,O_WAT_IN_C_id_lowerlayer[m])[1]/dt);
            }*/
            for (int f=frame;f<frame+1;f++){
                temp_O_X+=(nc_data.atom_coordinate(frame,O_WAT_IN_C_id_lowerlayer[m])[0]);
                //std::cout<<temp_O_X<<'\n'<<std::endl;
                temp_O_Y+=(nc_data.atom_coordinate(frame,O_WAT_IN_C_id_lowerlayer[m])[1]);
            }
            for (int i=0; i<x_points;i++){
                for(int j=0;j<y_points;j++){

                    density_matrix[i][j]+=density(X_DOWN+i*dx,Y_DOWN+j*dx,temp_O_X,temp_O_Y);//计算密度分布函数
                }
            }

        }
        for (int i=0; i<x_points;i++){
            for(int j=0;j<y_points;j++){
                if (density_matrix[i][j]<cavity_cut){
                    judge_matrix[i][j]=1;
                }
                else {
                    judge_matrix[i][j] = 0;
                }
                //通过密度矩阵判断是否为空穴
                //std::cout<<density_matrix[i][j]<<' ';
            }
            //std::cout<<std::endl;
        }
        int cavity_size=0;
        /*for (int i=0; i<x_points;i++){//挖掉边界
            cavity_size=0;
            calc_cavity_size(1,i,0,cavity_size,judge_matrix,calc_matrix,temp_coordinates,x_points,y_points);
            calc_cavity_size(1,i,y_points-1,cavity_size,judge_matrix,calc_matrix,temp_coordinates,x_points,y_points);
        }
        for (int j=0; j<y_points;j++){
            cavity_size=0;
            calc_cavity_size(1,0,j,cavity_size,judge_matrix,calc_matrix,temp_coordinates,x_points,y_points);
            calc_cavity_size(1,x_points-1,j,cavity_size,judge_matrix,calc_matrix,temp_coordinates,x_points,y_points);
        }*///别忘了改
        //std::cout<<'\n';
        std::cout<<frame_r<<' '<<std::endl;
        if (0<=x_r<x_points&&0<=y_r<y_points) {
            /*density_distribution[floor(density_matrix[x_r][y_r] / 0.0004)]++;
            std::cout << density_matrix[x_r][y_r] << std::endl;*/
            if (judge_matrix[x_r][y_r] == 1 && calc_matrix[x_r][y_r] == 0) {
                //target_jump++;
                int cavity_size=0;
                int cavity_size_1=0;
                temp_coordinates.clear();
                calc_cavity_size(1,x_r,y_r,cavity_size,judge_matrix,calc_matrix,temp_coordinates,x_points,y_points);
                //计算dt后的空穴变化
                frame+=dt;
                O_WAT_IN_C_id_lowerlayer.clear();
                O_WAT_IN_C_id_upperlayer.clear();
                for(index i = 0; i != O_WAT_id.size(); ++i)
                {
                    O_coor = nc_data.atom_coordinate(frame, O_WAT_id[i]);
                    if( O_coor[2] < Z_DOWN + (Z_UP - Z_DOWN)*7/15 && O_coor[2] >= Z_DOWN && O_coor[0] < X_UP && O_coor[0] > X_DOWN && O_coor[1] < Y_UP && O_coor[1] > Y_DOWN)//做修改，将水分为两层
                    {
                        O_WAT_IN_C_id_lowerlayer.push_back(O_WAT_id[i]);
                        //~ std::cout << O_WAT_id[i] << std::endl;
                    }
                    if( O_coor[2] <Z_UP && O_coor[2] >=  Z_UP - (Z_UP - Z_DOWN)*7/15  && O_coor[0] < X_UP && O_coor[0] > X_DOWN && O_coor[1] < Y_UP && O_coor[1] > Y_DOWN)//做修改，将水分为两层
                    {
                        O_WAT_IN_C_id_upperlayer.push_back(O_WAT_id[i]);
                        //~ std::cout << O_WAT_id[i] << std::endl;
                    }
                }
                for (int i = 0; i < x_points; i++) {
                    for (int j = 0; j < y_points; j++) {
                        density_matrix[i][j] = 0;
                        judge_matrix[i][j] = 0;
                        calc_matrix[i][j]=0;
                        //初始化
                    }
                }
                for(index m =0; m  != O_WAT_IN_C_id_lowerlayer.size(); ++m) {
                    temp_O_X=0;
                    temp_O_Y=0;
                    /*for (int f=frame;f<frame+dt;f++){
                        temp_O_X+=(nc_data.atom_coordinate(frame,O_WAT_IN_C_id_lowerlayer[m])[0]/dt);
                        //std::cout<<temp_O_X<<'\n'<<std::endl;
                        temp_O_Y+=(nc_data.atom_coordinate(frame,O_WAT_IN_C_id_lowerlayer[m])[1]/dt);
                    }*/
                    for (int f=frame;f<frame+1;f++){
                        temp_O_X+=(nc_data.atom_coordinate(frame,O_WAT_IN_C_id_lowerlayer[m])[0]);
                        //std::cout<<temp_O_X<<'\n'<<std::endl;
                        temp_O_Y+=(nc_data.atom_coordinate(frame,O_WAT_IN_C_id_lowerlayer[m])[1]);
                    }
                    for (int i=0; i<x_points;i++){
                        for(int j=0;j<y_points;j++){
                            density_matrix[i][j]+=density(X_DOWN+i*dx,Y_DOWN+j*dx,temp_O_X,temp_O_Y);//计算密度分布函数
                        }


                    }

                }
                for (int i=0; i<x_points;i++){
                    for(int j=0;j<y_points;j++){
                        if (density_matrix[i][j]<cavity_cut){
                            judge_matrix[i][j]=1;
                        }
                        else {
                            judge_matrix[i][j] = 0;
                        }
                        //通过密度矩阵判断是否为空穴
                        //std::cout<<judge_matrix[i][j]<<' ';
                    }
                    //std::cout<<std::endl;
                }
                for (;!temp_coordinates.empty();){
                    int *temp_xy = temp_coordinates.back();
                    temp_coordinates.pop_back();
                    if (judge_matrix[temp_xy[0]][temp_xy[1]] == 1) {
                        int temp_x = temp_xy[0];
                        int temp_y = temp_xy[1];
                        calc_cavity_size(1,temp_x,temp_y,cavity_size_1,judge_matrix,calc_matrix,temp_coordinates,x_points,y_points);

                        break;
                        std::cout<<temp_x<<' '<<temp_y<<' ';
                    }

                }
                double cavity_size_change=(cavity_size_1-cavity_size)*pow(dx,2);
                outfile<<frame_r<<std::setw(7)<<cavity_size*pow(dx,2)<<
                         std::setw(7)<<cavity_size_1*pow(dx,2)<<std::setw(7)<<cavity_size_change<<std::endl;
                size_change_distribution[floor(cavity_size_change / pow(dx,2))]++;
            } else {
                nt_jump++;
                /*std::cout<<nc<<' '<<frame<<' '<<x_r<<' '<<y_r<<std::endl;
                for (index m=0;m<O_WAT_IN_C_id_lowerlayer.size();m++){
                    std::cout<<O_WAT_IN_C_id_lowerlayer[m]<<' '<<O_WAT_IN_C_id_lowerlayer[m]+1<<' '<<O_WAT_IN_C_id_lowerlayer[m]+2<<' ';
                    std::cout<<O_WAT_IN_C_id_upperlayer[m]<<' '<<O_WAT_IN_C_id_upperlayer[m]+1<<' '<<O_WAT_IN_C_id_upperlayer[m]+2<<' ';
                }
                std::cout<<std::endl;*/
            }
        }
    }
    std::cout<<target_jump<<' '<<nt_jump<<std::endl;
    for (int i=0;i<3000;i++) {
        outfile2<<i*0.04<<std::setw(7)<<size_change_distribution[i]<<std::endl;
    }
    infile.close();
    outfile.close();
    outfile2.close();
    //infile1.close();
    return 0;
}
