
#include "amber_netcdf.hpp"
#include "amber_parm_1.0.hpp"
#include "vector_calc.h"
#define start_nc 37
#define end_nc  100
#define dt 4
#define frame_extend 100
#define name_parm7 "density_dis9a5.parm7"
//~ #define name_nc "water_ion_graphene_10a5"
typedef std::vector<double>::size_type index;

int main() {
    std::ifstream infile;
    std::ofstream outfile0;
    outfile0.open("jump/cutoff_1a3/transition_path_stat");
    static std::vector<double> cavity_frames;
    static std::vector<int> cavity_frame;
    char name_nc[64];
    sprintf(name_nc, "nc/density_dis9a5_%d.nc", start_nc);
    amber_parm parm_name(name_parm7);
    nctraj nc_data(name_nc);

    infile.open("jump/cutoff_1a3/transition_path_index_start_finish_down");
    while (!infile.fail()) {
        double num[4];
        infile >> num[0] >> num[1]>>num[2]>>num[3];
        int frame_r_st = num[1]-frame_extend;
        int frame_r_ed = num[2]+frame_extend;
        int Target_ID = num[0];
        for (int frame_r=frame_r_st;frame_r<frame_r_ed&&frame_r<1000000;frame_r+=dt) {
            int nc = start_nc;
            int total_frame = start_nc * 10000;
            while (total_frame <= frame_r) {
                sprintf(name_nc, "nc/density_dis9a5_%d.nc", nc);
                nctraj nc_data(name_nc);
                total_frame += nc_data.frames_number();
                nc++;
            }
            nc--;
            index frame = frame_r - total_frame + nc_data.frames_number();
            outfile0<<Target_ID<<std::setw(10)<<nc<<std::setw(10)<<frame<<std::setw(10)<<num[3]<<std::endl;
        }
    }
    infile.close();
    outfile0.close();
    //infile1.close();
    return 0;
}
