//
// Created by utena on 18-3-1.
//
//
// Created by utena on 17-10-27.
//
#include "Environment& Algorithm.hpp"
#include <string>
#define max_transition_time 20000
#define tot_frame 20000
#define frame_segment 400
#define max_cluster_size 24
#define deviation_y_center 30
#define deviation_x_center 30
#define start_nc 0
#define end_nc  0
#define z_axis_num 200
#define frame_interval 0.25
using std::cout;
using std::vector;
using std::ifstream;
using std::cin;
using std::string;
using std::endl;
using std::setw;

int main(int argc,char* argv[]) {
    CreatDir("z_axis_diffusion");
    std::cout << "program to calculate diffusion along z axis" << "\n" << std::endl;
    //cout<<"Loading trajectories..."<<endl;

    std::vector<double> O_coor, O1_coor,O2_coor, H1O_coor, H2O_coor,H1_coor, H2_coor,O_coor_0;
    double Diffusion_sum_z[frame_segment]={0};
    double Diffusion_sum_xy[frame_segment]={0};
    int WAT_ID_to_i[12000];
    for (int i = 0; i < 12000; i++) WAT_ID_to_i[i] = -1;
    char dist[16];
    char frcmod[16];
    sprintf(dist,"%s",argv[1]);
    sprintf(frcmod,"%s",argv[2]);
    int heat_temp=std::stoi(argv[3]);
    int equ_temp=std::stoi(argv[4]);
    typedef std::vector<double>::size_type index;
    char file_prefix[128];
    sprintf(file_prefix,"../dis%s_%s_WATBOX/MD_HEAT%dK_EQ%dK/Equilibrium/density_dis%s_WATBOX",dist,frcmod,heat_temp,equ_temp,dist);
    char name_parm7[128];
    sprintf(name_parm7, "%s.parm7", file_prefix);
    std::ofstream outfile;
    char name_output[128];
    sprintf(name_output, "z_axis_diffusion/%s.dat", dist);
    outfile.open(name_output);
    std::cout << "calculate XYZ_limit" << std::endl;
    std::vector<double> z_axis_count(z_axis_num, 0);
    double C_z_coor_sum_1 = 0, C_z_coor_sum_2 = 0;
    double C_z_coor_average_1, C_z_coor_average_2;
    double C_x_coor_sum = 0, C_y_coor_sum = 0;
    double C_x_coor_center, C_y_coor_center;
    double Z_UP, Z_DOWN, Y_UP, Y_DOWN, X_UP, X_DOWN;
    amber_parm parm_nam(name_parm7);
    char name_nc[128];
    sprintf(name_nc, "%s_NVT.nc", file_prefix);
    nctraj data_nc(name_nc);
    for (index C_index = 0; C_index != C_NUM; ++C_index) {
        C_z_coor_sum_1 += data_nc.atom_coordinate(0, C_index)[2];
        C_z_coor_sum_2 += data_nc.atom_coordinate(0, (C_NUM + C_index))[2];
        C_y_coor_sum += data_nc.atom_coordinate(0, C_index)[1];
        C_x_coor_sum += data_nc.atom_coordinate(0, C_index)[0];
    }
    C_z_coor_average_1 = C_z_coor_sum_1/(C_NUM*(end_nc-start_nc+1));
    C_z_coor_average_2 = C_z_coor_sum_2/(C_NUM*(end_nc-start_nc+1));
    C_y_coor_center = C_y_coor_sum/(C_NUM*(end_nc-start_nc+1));
    C_x_coor_center = C_x_coor_sum/(C_NUM*(end_nc-start_nc+1));

    Z_UP = C_z_coor_average_2;
    Z_DOWN = C_z_coor_average_1;
    Y_UP = C_y_coor_center + deviation_y_center;
    Y_DOWN = C_y_coor_center - deviation_y_center;
    X_UP = C_x_coor_center + deviation_x_center;
    X_DOWN = C_x_coor_center - deviation_x_center;
    double dz=(Z_UP-Z_DOWN)/z_axis_num;

    std::cout << "Z_UP: " << Z_UP << std::endl;
    std::cout << "Z_DOWN: " << Z_DOWN << std::endl;
    std::cout << "X_UP: " << X_UP << std::endl;
    std::cout << "X_DOWN: " << X_DOWN << std::endl;
    std::cout << "Y_UP: " << Y_UP << std::endl;
    std::cout << "Y_DOWN: " << Y_DOWN << std::endl;
    std::cout << "dz: " << dz << std::endl;

    amber_parm parm_name(name_parm7);
    sprintf(name_nc, "%s_NVT.nc", file_prefix);
    nctraj nc_data(name_nc);
    int total_frame = nc_data.frames_number();
    std::vector<index> O_WAT_id = parm_name.id_by_type("OW");
    std::vector<index> O_WAT_IN_BOX_id;
    int tot_wat_num=O_WAT_id.size();
    int frame = 0;
    int WATcount=0;
    std::vector<double> o_coor;

    while (frame < total_frame-frame_segment) {
        if (frame%25 == 0) {
            O_WAT_IN_BOX_id.clear();
            for (index i = 0; i != O_WAT_id.size(); ++i) {
                o_coor = nc_data.atom_coordinate(frame, O_WAT_id[i]);
                if (o_coor[2] < Z_UP && o_coor[2] > Z_DOWN && o_coor[0] < X_UP && o_coor[0] > X_DOWN &&
                    o_coor[1] < Y_UP && o_coor[1] > Y_DOWN) {
                    O_WAT_IN_BOX_id.push_back(O_WAT_id[i]);
                    WATcount += 1;
                }
            }
            for (index i=0; i!=O_WAT_IN_BOX_id.size();i++){
                O_coor_0=nc_data.atom_coordinate(frame,O_WAT_IN_BOX_id[i]);
                for (int frame_diff=0;frame_diff<frame_segment;frame_diff++){
                    O_coor=nc_data.atom_coordinate(frame+frame_diff,O_WAT_IN_BOX_id[i]);
                    Diffusion_sum_z[frame_diff]+=pow((O_coor[2]-O_coor_0[2]),2);
                    Diffusion_sum_xy[frame_diff]+=(pow((O_coor[1]-O_coor_0[1]),2),pow((O_coor[0]-O_coor_0[0]),2));

                }
            }
            for (int i = 0; i != frame_segment; ++i) {
                //std::cout << i * dt << std::setw(15) << Diffusion_sum_xy[i] / WATcount << std::setw(15)
                          //<< Diffusion_sum_z[i] / WATcount << std::endl;
            }
        }
        frame++;
    }
    for (int i = 0; i != frame_segment; ++i) {
        outfile <<  i*frame_interval << std::setw(15) << Diffusion_sum_xy[i]/WATcount << std::setw(15) << Diffusion_sum_z[i]/WATcount << std::endl;
    }
    outfile.close();

    /*for (int WAT_ID_ref = 0; WAT_ID_ref != tot_wat_num; ++WAT_ID_ref) {
        char outfile_name[64];
        std::ofstream tracing_out;
        sprintf(outfile_name, "Analysis/Cluster_size/WAT_%d_nc_%02d-%02d.dat",O_WAT_id[WAT_ID_ref], start_nc,end_nc);
        tracing_out.open(outfile_name);
        for (int t=0;t!=Cluster_Size_by_WATID_frame[0].size();t++){
            tracing_out<<t*0.25<<setw(10)<<Cluster_Size_by_WATID_frame[WAT_ID_ref][t]<<endl;
        }
        tracing_out.close();
    }*/
    cout<<"Loading complete"<<endl;


    return 0;
}

