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
#define deviation_y_center 10
#define deviation_x_center 10
#define start_nc 49
#define end_nc  100
#define dt 1
#define relaxation_time 25
#define name_parm7 "density_dis9a5.parm7"
#define psi 2.4
#define cavity_cut 0.006
#define dx 0.2
#define cut_off 24
#define max_size 2000
#define max_life 10000
#define max_circumference 400
#define rapid_change_constant 5

double density(double x,double y, double x_W, double y_W){
    double r= sqrt(pow((x-x_W),2)+pow((y-y_W),2));
    return (1/(2*M_PI*pow(psi,2)))*exp(-2*(pow((r/psi),2)));
}

int calc_cavity_size(int &circumference,int x,int y,int &last_step_size, int **judge_matrix, int **calc_matrix,std::vector<int*> &temp_cavity_coordinate,int x_points, int y_points){
    int *coordinate_vector=new int [2];
	if (calc_matrix[x][y]==0){
        if (x==0 ||x==x_points-1||y==0||y==y_points-1){
            circumference++;
        } else{
            if (!(judge_matrix[x-1][y]==1&judge_matrix[x+1][y]==1&judge_matrix[x][y-1]==1&judge_matrix[x][y+1]==1)){
                circumference++;
            }
        }
        last_step_size+=1;
        calc_matrix[x][y]=1;
        if (x-1>=0){
            if (judge_matrix[x-1][y]==1 && calc_matrix[x-1][y]==0){
				coordinate_vector[0]=x;
				coordinate_vector[1]=y;
				temp_cavity_coordinate.push_back(coordinate_vector);
                calc_cavity_size(circumference,x-1,y,last_step_size,judge_matrix,calc_matrix,temp_cavity_coordinate,x_points,y_points);
            }
        }
        if (y-1>=0){
            if (judge_matrix[x][y-1]==1 && calc_matrix[x][y-1]==0){
				coordinate_vector[0]=x;
				coordinate_vector[1]=y;
				temp_cavity_coordinate.push_back(coordinate_vector);
                calc_cavity_size(circumference,x,y-1,last_step_size,judge_matrix,calc_matrix,temp_cavity_coordinate,x_points,y_points);
            }
        }
        if (x+1<x_points){
            if (judge_matrix[x+1][y]==1 && calc_matrix[x+1][y]==0){
				coordinate_vector[0]=x;
				coordinate_vector[1]=y;
				temp_cavity_coordinate.push_back(coordinate_vector);
                calc_cavity_size(circumference,x+1,y,last_step_size,judge_matrix,calc_matrix,temp_cavity_coordinate,x_points,y_points);
            }
        }
        if (y+1<y_points){
            if (judge_matrix[x][y+1]==1 && calc_matrix[x][y+1]==0){
				coordinate_vector[0]=x;
				coordinate_vector[1]=y;
				temp_cavity_coordinate.push_back(coordinate_vector);
                calc_cavity_size(circumference,x,y+1,last_step_size,judge_matrix,calc_matrix,temp_cavity_coordinate,x_points,y_points);
            }
        }
    }
    return last_step_size;
}
int main() {
    std::ofstream outfile;
    std::ofstream outfile2;
    outfile.open("Rapid_size_change");

    typedef std::vector<double>::size_type index;
    std::cout << "program to find the cavities" << "\n" << std::endl;
    int x_points = deviation_x_center * 2 / dx;
    int y_points = deviation_y_center * 2 / dx;
	double max_cavity_index=0;
    std::vector<double> cavity_quantities_by_size(max_size,0);
    std::vector<double> cavity_quantities_by_circumference(max_circumference,0);
    std::vector<double> cavity_quantities_by_life(max_life,0);

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
    for (int nc = start_nc; nc != end_nc + 1; ++nc) {
        char name_nc[64];
        char name_out2[64];
        sprintf(name_nc, "density_dis9a5_%d.nc", nc);
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


        double C_z_coor_sum_1 = 0, C_z_coor_sum_2 = 0;
        double C_z_coor_average_1, C_z_coor_average_2;
        double C_x_coor_sum = 0, C_y_coor_sum = 0;
        double C_x_coor_center, C_y_coor_center;

        double Z_UP, Z_DOWN, Y_UP, Y_DOWN, X_UP, X_DOWN;

        for (int frame = 0; frame < total_frame; frame+=relaxation_time)
            //~ for(int frame =0; frame != 2       ;++frame)
        {
            if ((frame % 1000) == 0) {
                std::cout << "frame_now : " << frame << std::endl;
            }
            if (frame == 0) {
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
                for (int i=0; i<x_points;i++) {
                    for (int j = 0; j < y_points; j++) {
                        sum_matrix[i][j] = 0;
                    }
                }
            }
            std::vector<double> O_coor,O_coor_initial;
            if(frame %dt==0 && frame!=total_frame) {
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
                for (int frame_t=frame;frame_t<frame+relaxation_time;frame_t+=dt) {
                    for (index i = 0; i != O_WAT_id.size(); ++i) {
                        O_coor = nc_data.atom_coordinate(frame, O_WAT_id[i]);
                        if (O_coor[2] < Z_DOWN + (Z_UP - Z_DOWN) * 7 / 15 && O_coor[2] > Z_DOWN && O_coor[0] < X_UP &&
                            O_coor[0] > X_DOWN && O_coor[1] < Y_UP && O_coor[1] > Y_DOWN)//做修改，将水分为两层
                        {
                            O_WAT_IN_C_id_lowerlayer.push_back(O_WAT_id[i]);
                            //~ std::cout << O_WAT_id[i] << std::endl;
                        }
                        if (O_coor[2] < Z_UP && O_coor[2] > Z_UP - (Z_UP - Z_DOWN) * 7 / 15 && O_coor[0] < X_UP &&
                            O_coor[0] > X_DOWN && O_coor[1] < Y_UP && O_coor[1] > Y_DOWN)//做修改，将水分为两层
                        {
                            O_WAT_IN_C_id_upperlayer.push_back(O_WAT_id[i]);
                            //~ std::cout << O_WAT_id[i] << std::endl;
                        }
                    }
                    double temp_O_X;
                    double temp_O_Y;
                    for (index m = 0; m != O_WAT_IN_C_id_lowerlayer.size(); ++m) {
                        temp_O_X = 0;
                        temp_O_Y = 0;
                        for (int f = frame_t; f < frame_t + dt; f++) {
                            temp_O_X += (nc_data.atom_coordinate(f, O_WAT_IN_C_id_lowerlayer[m])[0] / dt);
                            //std::cout<<temp_O_X<<'\n'<<std::endl;
                            temp_O_Y += (nc_data.atom_coordinate(f, O_WAT_IN_C_id_lowerlayer[m])[1] / dt);
                        }
                        for (int i = 0; i < x_points; i++) {
                            for (int j = 0; j < y_points; j++) {
                                density_matrix[i][j] += density(X_DOWN + i * dx, Y_DOWN + j * dx, temp_O_X,
                                                                temp_O_Y)*dt/relaxation_time;//计算密度分布函数
                            }


                        }
                    }
                }
                std::cout<<"step: "<<frame/dt<<'\n'<<std::endl;
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

                int large_cavity_num=0;
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
				temp_coordinates.clear();
                bool is_active= false;
                new_active_index.clear();
                //std::cout<<last_step_cavity_coordinate_by_index.size()<<std::endl;
                /*for (int i=0;i<active_cavity_index.size();i++){//计算空穴覆盖情况
                    is_active=false;

                    //std::cout<<active_cavity_index[i]<<' '<<last_step_cavity_coordinate_by_index[active_cavity_index[i]].size()<<std::endl;

                    int temp_x=0;int temp_y=0;
                    for (;!last_step_cavity_coordinate_by_index[active_cavity_index[i]].empty();){
                        int *temp_xy = last_step_cavity_coordinate_by_index[active_cavity_index[i]].back();
                        last_step_cavity_coordinate_by_index[active_cavity_index[i]].pop_back();
                        if (judge_matrix[temp_xy[0]][temp_xy[1]] == 1) {
                                is_active = true;//版本1 有重合就是同一个
                                temp_x = temp_xy[0];
                                temp_y = temp_xy[1];
                                //std::cout<<temp_x<<' '<<temp_y<<' ';
                        }

                    }
                    last_step_cavity_coordinate_by_index[active_cavity_index[i]].clear();
                    cavity_size=0;
                    if (is_active){

                        calc_cavity_size(active_cavity_index[i],temp_x,temp_y,cavity_size,judge_matrix,calc_matrix,last_step_cavity_coordinate_by_index[active_cavity_index[i]],x_points,y_points);

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
                 */
                for (int i=0; i<x_points;i++){//深搜统计新的空穴大小和数目
                    for(int j=0;j<y_points;j++){
                        if (judge_matrix[i][j]==1&& calc_matrix[i][j]==0){
                            cavity_size=0;
                            cavity_circ=0;
                            temp_coordinates.clear();
                            calc_cavity_size(cavity_circ,i,j,cavity_size,judge_matrix,calc_matrix,temp_coordinates,x_points,y_points);//
                            //temp_coordinates.clear();//完全忽略面积信息，内存占用极大时使用
                            //std::cout<<temp_coordinates.size()<<std::endl;
                            //last_step_cavity_coordinate_by_index.push_back(temp_coordinates);
                            std::vector<int> init_cavity_size_vector;
                            //init_cavity_size_vector.push_back(cavity_size);
                            //cavity_size_by_time.push_back(init_cavity_size_vector);
                            //cavity_frame.push_back(frame);
                            //active_cavity_index.push_back(cavity_num);
                            cavity_quantities_by_size[cavity_size]++;
                            cavity_quantities_by_circumference[cavity_circ]++;
                            cavity_num++;
                            if (cavity_size*dx*dx>=cut_off){
                                large_cavity_num++;
                                /*for (index m=0;m<O_WAT_IN_C_id_lowerlayer.size();m++){
                                    //std::cout<<O_WAT_IN_C_id_lowerlayer[m]<<' '<<O_WAT_IN_C_id_lowerlayer[m]+1<<' '<<O_WAT_IN_C_id_lowerlayer[m]+2<<' ';
                                }*/
                                //std::cout<<std::endl;
                            }
                            //std::cout<<cavity_size<<std::endl;
                        }
                    }
                }
                //std::cout<<frame <<' '<<large_cavity_num<<std::endl;

            }
        }
        std::ofstream outfile1;
        char out_name[64];
        sprintf(out_name,"cavity_size_density_relaxation_%dfs.dat",relaxation_time*4);
        outfile1.open(out_name);
        for (index i=0;i<cavity_quantities_by_size.size();i++){
            std::cout<<i*pow(dx,2)<<std::setw(15)<<cavity_quantities_by_size[i]/((nc-start_nc+1)*10000/relaxation_time)<<std::endl;
            outfile1<<i*pow(dx,2)<<std::setw(15)<<cavity_quantities_by_size[i]/((nc-start_nc+1)*10000/relaxation_time)<<std::endl;
        }
        for (index i=0; i<max_life;i++) {

        }
        outfile1.close();
        std::ofstream outfile3;
        char out3_name[64];
        sprintf(out3_name,"cavity_circumference_density_relaxation_%dfs.dat",relaxation_time*4);
        outfile3.open(out3_name);
        for (index i=0;i<max_circumference;i++){
            std::cout<<i*dx<<std::setw(15)<<cavity_quantities_by_circumference[i]/((nc-start_nc+1)*10000/relaxation_time)<<std::endl;
            outfile3<<i*dx<<std::setw(15)<<cavity_quantities_by_circumference[i]/((nc-start_nc+1)*10000/relaxation_time)<<std::endl;
        }
        for (index i=0; i<max_life;i++) {

        }
        outfile3.close();
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


    outfile.close();
    return 0;
}
