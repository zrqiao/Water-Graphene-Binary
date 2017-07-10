
//
// Created by utena on 17-5-8.
//
#include "enthaply/amber_netcdf.hpp"
#include "enthaply/amber_parm_1.0.hpp"
#include<cmath>
#include<cstdlib>
#include<cstring>
#include<cfloat>
#define z_axis_modify  2
#define deviation_y_center 26
#define deviation_x_center 26
#define start_nc 37
#define end_nc  52
#define dt 4
#define relaxation_time 4
#define name_parm7 "nc/density_dis9a5.parm7"
#define psi 2.4
#define cavity_cut 0.006
#define dx 0.65
#define cut_off 24
#define max_size 2000
#define max_life 10000
#define max_circumference 400
#define z_points    160
#define large_cutoff 5.4
#define rapid_change_constant 5
static int x_points = deviation_x_center * 2 / dx;
static int y_points = deviation_y_center * 2 / dx;

double density(double x,double y, double x_W, double y_W){
    double r= sqrt(pow((x-x_W),2)+pow((y-y_W),2));
    return (1/(2*M_PI*pow(psi,2)))*exp(-2*(pow((r/psi),2)));
}

int calc_cavity_size(int &last_step_size, int &circumference, int x, int y, int **judge_matrix, int **calc_matrix) {
    int *coordinate_vector=new int [2];
        if (judge_matrix [x][y]==1 && calc_matrix[x][y] == 0) {
            if (x == 0 || x == x_points - 1 || y == 0 || y == y_points - 1) {
                circumference++;
            } else {
                if (!(judge_matrix[x - 1][y] == 1 & judge_matrix[x + 1][y] == 1 & judge_matrix[x][y - 1] == 1 &
                      judge_matrix[x][y + 1] == 1)) {
                    circumference++;
                }
            }
            last_step_size += 1;
            calc_matrix[x][y] = 1;
            if (x - 1 >= 0) {
                if (judge_matrix[x - 1][y] == 1 && calc_matrix[x - 1][y] == 0) {
                    coordinate_vector[0] = x;
                    coordinate_vector[1] = y;
                    calc_cavity_size(last_step_size, circumference, x - 1, y, judge_matrix, calc_matrix);
                }
            }
            if (y - 1 >= 0) {
                if (judge_matrix[x][y - 1] == 1 && calc_matrix[x][y - 1] == 0) {
                    coordinate_vector[0] = x;
                    coordinate_vector[1] = y;
                    calc_cavity_size(last_step_size, circumference, x, y - 1, judge_matrix, calc_matrix);
                }
            }
            if (x + 1 < x_points) {
                if (judge_matrix[x + 1][y] == 1 && calc_matrix[x + 1][y] == 0) {
                    coordinate_vector[0] = x;
                    coordinate_vector[1] = y;
                    calc_cavity_size(last_step_size, circumference, x + 1, y, judge_matrix, calc_matrix);
                }
            }
            if (y + 1 < y_points) {
                if (judge_matrix[x][y + 1] == 1 && calc_matrix[x][y + 1] == 0) {
                    coordinate_vector[0] = x;
                    coordinate_vector[1] = y;
                    calc_cavity_size(last_step_size, circumference, x, y + 1, judge_matrix, calc_matrix);
                }
            }
        }
    return last_step_size;
}
int main() {
    std::ifstream infile;
    infile.open("jump/cutoff_1a3/transiton_stat_fastinput");
    typedef std::vector<double>::size_type index;
    std::cout << "program to find size-z correlation" << "\n" << std::endl;

    double max_cavity_index=0;
    std::vector<double> cavity_size_sum(z_points,0);
    std::vector<double> large_cavity_probability(z_points,0);
    std::vector<double> total_wat_z(z_points,0);
    std::vector<int> active_cavity_index;
    std::vector<int> new_active_index;
    std::vector<std::vector<int*>> last_step_cavity_coordinate_by_index;
    static std::vector<std::vector<int>> cavity_size_by_time;
    static std::vector<int> cavity_frame;
    static int cavity_num=0;
    static double **density_matrix = new double *[x_points];
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
    static double **sum_matrix = new double *[x_points];
    for (int i = 0; i < x_points; i++) {
        sum_matrix[i] = new double[y_points];
    }
    static std::vector<double> cavity_by_index;
    std::vector<int*> temp_coordinates;
    double num[4];
    int Target_ID, direction;
    int frame_read, nc_read;
    while (num[1]!=start_nc) {
        infile >> num[0] >> num[1] >> num[2] >> num[3];
        Target_ID = num[0];
        nc_read = num[1];
        frame_read = num[2];
        direction = num[3];
    }
    for (int nc = start_nc; nc != end_nc + 1; ++nc) {
        char name_nc[64];
        char name_out2[64];
        sprintf(name_nc, "nc/density_dis9a5_%d.nc", nc);
        amber_parm parm_name(name_parm7);
        nctraj nc_data(name_nc);
        printf("nc_file = %d\n", nc);
        int total_frame = nc_data.frames_number();
        std::cout << "total frame: " << total_frame << std::endl;
        std::vector<index> O_WAT_id = parm_name.id_by_type("OW");
        std::cout << "total water:  " << O_WAT_id.size() << "\n" << std::endl;
        std::vector<index> O_WAT_IN_C_id_upperlayer;
        std::vector<index> O_WAT_IN_C_id_lowerlayer;

        double dz;

        double Z_UP, Z_DOWN, Y_UP, Y_DOWN, X_UP, X_DOWN;
        std::vector<double> O1_coor, O2_coor, H1_coor, H2_coor;

        for (int frame = 0; frame < total_frame; frame+=relaxation_time){
            if ((frame % 100) == 0) {
                double C_z_coor_sum_1 = 0, C_z_coor_sum_2 = 0;
                double C_z_coor_average_1, C_z_coor_average_2;
                double C_x_coor_sum = 0, C_y_coor_sum = 0;
                double C_x_coor_center, C_y_coor_center;
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
                dz=(Z_UP-Z_DOWN)/z_points;
            }
            std::vector<double> O_coor,O_coor_initial;
            if(frame!=total_frame && ( nc_read==nc && frame_read==frame)) {
                O_WAT_IN_C_id_lowerlayer.clear();
                O_WAT_IN_C_id_upperlayer.clear();
                for (int i = 0; i < x_points; i++) {
                    for (int j = 0; j < y_points; j++) {
                        density_matrix[i][j] = 0;
                        judge_matrix[i][j] = 0;
                        calc_matrix[i][j] = 0;
                        //初始化
                    }
                }
                for (index i = 0; i != O_WAT_id.size(); ++i) {
                    O_coor = nc_data.atom_coordinate(frame, O_WAT_id[i]);
                    if (O_coor[2] <  (Z_UP + Z_DOWN) /2 -0.8 && O_coor[2] > Z_DOWN && O_coor[0] < X_UP &&
                        O_coor[0] > X_DOWN && O_coor[1] < Y_UP && O_coor[1] > Y_DOWN)//做修改，将水分为两层
                    {
                        O_WAT_IN_C_id_lowerlayer.push_back(O_WAT_id[i]);
                        //~ std::cout << O_WAT_id[i] << std::endl;
                    }
                    if (O_coor[2] < Z_UP && O_coor[2] > 0.8+ (Z_UP + Z_DOWN) /2 && O_coor[0] < X_UP &&
                        O_coor[0] > X_DOWN && O_coor[1] < Y_UP && O_coor[1] > Y_DOWN)//做修改，将水分为两层
                    {
                        O_WAT_IN_C_id_upperlayer.push_back(O_WAT_id[i]);
                        //~ std::cout << O_WAT_id[i] << std::endl;
                    }
                }
                for (int frame_t=frame;frame_t<frame+relaxation_time;frame_t+=dt) {
                    double temp_O_X;
                    double temp_O_Y;
                    for (index m = 0; m != O_WAT_IN_C_id_upperlayer.size(); ++m) {
                        temp_O_X = 0;
                        temp_O_Y = 0;
                        for (int f = frame_t; f < frame_t + dt; f++) {
                            temp_O_X += (nc_data.atom_coordinate(f, O_WAT_IN_C_id_upperlayer[m])[0] / dt);
                            //std::cout<<temp_O_X<<'\n'<<std::endl;
                            temp_O_Y += (nc_data.atom_coordinate(f, O_WAT_IN_C_id_upperlayer[m])[1] / dt);
                        }
                        for (int i = 0; i < x_points; i++) {
                            for (int j = 0; j < y_points; j++) {
                                density_matrix[i][j] += density(X_DOWN + i * dx, Y_DOWN + j * dx, temp_O_X,
                                                                temp_O_Y)*dt/relaxation_time;//计算密度分布函数
                            }
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
                        //sum_matrix[i][j]+=judge_matrix[i][j];
                    }
                    //std::cout<<std::endl;
                }
                int cavity_size, cavity_circ;
                /*for (int i=0; i<x_points;i++){//挖掉边界
                    cavity_size=0;
                    calc_cavity_size(1,i,0,cavity_size,judge_matrix,calc_matrix,temp_coordinates,x_points,y_points);
                    calc_cavity_size(1,i,y_points-1,cavity_size,judge_matrix,calc_matrix,temp_coordinates,x_points,y_points);
                }
                for (int j=0; j<y_points;j++){
                    cavity_size=0;
                    calc_cavity_size(1,0,j,cavity_size,judge_matrix,calc_matrix,temp_coordinates,x_points,y_points);
                    calc_cavity_size(1,x_points-1,j,cavity_size,judge_matrix,calc_matrix,temp_coordinates,x_points,y_points);
                }*/
                while (nc_read==nc && frame_read==frame) {
                    std::cout<<Target_ID<<std::setw(10)<<nc_read<<std::setw(10)<<frame_read<<std::setw(10)<<direction<<std::endl;
                    if (direction==1) {//只考虑自下向上
                        for (int i = 0; i < x_points; i++) {
                            for (int j = 0; j < y_points; j++) {
                                calc_matrix[i][j] = 0;
                            }
                        }
                        O1_coor = nc_data.atom_coordinate(frame, Target_ID);
                        int x=int(round((O1_coor[0] - X_DOWN) / dx)),y=int(round((O1_coor[1] - Y_DOWN) / dx));
                        if (0<=x&&x<x_points && 0<=y&&y<y_points) {
                            //std::cout << x <<std::setw(5)<<y<< std::endl;
                            cavity_circ = 0;
                            cavity_size = 0;
                            calc_cavity_size(cavity_circ, cavity_size, x,y , judge_matrix, calc_matrix);
                            double size = cavity_size * pow(dx, 2);
                            total_wat_z[int(round((O1_coor[2] - Z_DOWN) / (dz)))]++;
                            cavity_size_sum[int(round((O1_coor[2] - Z_DOWN) / (dz)))] += size;
                            std::cout<<size<<std::endl;
                            if (size >= large_cutoff)
                                large_cavity_probability[int(round((O1_coor[2] - Z_DOWN) / (dz)))]++;
                        }
                    }
                    if (!infile.fail()){
                        double num[4];
                        infile >> num[0] >> num[1]>>num[2]>>num[3];
                        Target_ID=num[0];
                        nc_read=num[1];
                        frame_read=num[2];
                        direction=num[3];

                    }
                }
                //std::cout<<frame <<' '<<large_cavity_num<<std::endl;

            }
        }
        std::ofstream outfile;
        std::ofstream outfile2;
        outfile.open("cavity_jump_correlation/cutoff_1a3_average_cavitysize_by_z");
        outfile2.open("cavity_jump_correlation/cutoff_1a3_largecavity_prob_by_z");
        for(int i =1; i !=z_points; ++i)
        {
            outfile<<Z_DOWN+i*dz<<std::setw(14)<<cavity_size_sum[i]/total_wat_z[i]<<std::endl;
            outfile2<<Z_DOWN+i*dz<<std::setw(14)<<large_cavity_probability[i]/total_wat_z[i]<<std::endl;
        }
        outfile.close();
        outfile2.close();
    }
    /*outfile.open("calc_cavity");
    for (index i=0; i<max_size;i++){
        std::cout<<i*pow(dx,2)<<std::setw(10)<<cavity_quantities_by_size[i]<<std::endl;
        outfile<<i*pow(dx,2)<<std::setw(10)<<cavity_quantities_by_size[i]<<std::endl;
    }*/
    /*for (int i=0; i<x_points;i++) {
        for (int j = 0; j < y_points; j++) {
            std::cout << i * dx << std::setw(10) << j * dx << std::setw(10) << sum_matrix[i][j] << std::endl;
        }
    }*/
    return 0;
}
