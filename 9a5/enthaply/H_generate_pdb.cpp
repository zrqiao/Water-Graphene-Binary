#include "amber_netcdf.hpp"
#include "amber_parm_1.0.hpp"
#include "vector_calc.h"


#define z_axis_modify  2
#define deviation_y_center 20
#define deviation_x_center 20
#define frame 1000
#define num_wat 100
#define name_parm7 "density_dis6a5.parm7" 
#define name_nc "density_dis6a5_15.nc"
//~ #define name_nc "water_ion_graphene_10a5"

int main()
{
	typedef std::vector<double>::size_type index;
	std::cout << "program to find the 100 water index"<< "\n"<< std::endl;
	std::cout << "calculate XYZ_limit" << std::endl;
	
    double C_z_coor_sum_1 = 0, C_z_coor_sum_2 = 0;
    double C_z_coor_average_1, C_z_coor_average_2;	
	double C_x_coor_sum = 0, C_y_coor_sum = 0;
    double C_x_coor_center, C_y_coor_center;     
    double Z_UP,Z_DOWN,Y_UP,Y_DOWN,X_UP,X_DOWN;
    //~ X_UP = -1000;
    //~ Y_UP = -1000;
    //~ X_DOWN = 1000;
    //~ Y_DOWN = 1000;
    	
	amber_parm parm_name(name_parm7);
	nctraj nc_data(name_nc);
	
	std::cout << "cell length:  "<< std::setw(15)<<nc_data.cell_lengths(frame)[0]
	                                                    << std::setw(15)<<nc_data.cell_lengths(frame)[1]
	                                                    << std::setw(15)<<nc_data.cell_lengths(frame)[2]  << std::endl;
	
	

	for( index C_index = 0 ; C_index != 1400 ; ++C_index)
	{
		       C_z_coor_sum_1 += nc_data.atom_coordinate(0, C_index)[2];
		       C_z_coor_sum_2 += nc_data.atom_coordinate(0, (1400+C_index))[2];
		       C_y_coor_sum += nc_data.atom_coordinate(0, C_index)[1];
   		       C_x_coor_sum += nc_data.atom_coordinate(0, C_index)[0];	
   		       
   		       //~ std::cout << nc_data.atom_coordinate(0, C_index)[0]<< std::endl;
   		       
   		       
   		       //~ if(data_nc.atom_coordinate(0, C_index)[1] > Y_UP)Y_UP=data_nc.atom_coordinate(0, C_index)[1];
      		   //~ if(data_nc.atom_coordinate(0, C_index)[1] < Y_DOWN)Y_DOWN=data_nc.atom_coordinate(0, C_index)[1];		       
   		       //~ if(data_nc.atom_coordinate(0, C_index)[0] > X_UP)X_UP=data_nc.atom_coordinate(0, C_index)[0];   		
    		   //~ if(data_nc.atom_coordinate(0, C_index)[0] < X_DOWN)X_DOWN=data_nc.atom_coordinate(0, C_index)[0];   	  		              
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
    
    
    std::vector<std::vector<index>> O_H1_H2_WAT_IN_C_id;
    std::vector<index> O_WAT_IN_C_id, H1_WAT_IN_C_id, H2_WAT_IN_C_id;
    std::vector<std::vector<double>> O_WAT_IN_C_coor, H1_WAT_IN_C_coor, H2_WAT_IN_C_coor;
    std::vector<double> O_WAT_IN_C_dis;
    std::vector<double> O_coor, H1_coor, H2_coor;
    std::vector<index> O_WAT_id = parm_name.id_by_type("OW");
    std::vector<index> graphene_id = parm_name.ids_by_resname("gra");
    
    //~ for(index i = 0; i != graphene_id.size(); ++i)
    //~ {
		//~ std::cout << graphene_id[i]<< std::endl;
	//~ }
       
    
	for( index i = 0 ; i != O_WAT_id.size() ; ++i)                                                             
	{
	    O_coor = nc_data.atom_coordinate(frame,O_WAT_id[i]);
		if( O_coor[2] < Z_UP && O_coor[2] > Z_DOWN && O_coor[0] < X_UP && O_coor[0] > X_DOWN && O_coor[1] < Y_UP && O_coor[1] > Y_DOWN)
		{
			//~ if((pow((O_coor[0]-C_x_coor_center),2)+pow((O_coor[1]-C_y_coor_center),2)) <  max_distance)
			//~ {
				//~ index_o_wat[largest_location] = O_WAT_id[i];
		    //~ }
		    double distance;
		    distance = pow(O_coor[0]-C_x_coor_center,2)+pow(O_coor[1]-C_y_coor_center,2);
		    O_WAT_IN_C_dis.push_back(distance);
		    O_WAT_IN_C_coor.push_back(O_coor);
		    H1_WAT_IN_C_coor.push_back(nc_data.atom_coordinate(frame,O_WAT_id[i]+1));
		    H2_WAT_IN_C_coor.push_back(nc_data.atom_coordinate(frame,O_WAT_id[i]+2));		    
			O_WAT_IN_C_id.push_back(O_WAT_id[i]);
			H1_WAT_IN_C_id.push_back(O_WAT_id[i]+1);
			H2_WAT_IN_C_id.push_back(O_WAT_id[i]+2);
		}
	 }
	 
	 O_H1_H2_WAT_IN_C_id.push_back(O_WAT_IN_C_id);
	 O_H1_H2_WAT_IN_C_id.push_back(H1_WAT_IN_C_id);
	 O_H1_H2_WAT_IN_C_id.push_back(H2_WAT_IN_C_id);
	 
	 double temp1;
	 index temp2;
	 std::vector<double> temp3;
	 
	 
	 for(index i = 0; i != (O_WAT_IN_C_id.size()-1); ++i)
	 {
		 for(index j = i+1; j != O_WAT_IN_C_id.size(); ++j)
		 {
			 if( O_WAT_IN_C_dis[i]> O_WAT_IN_C_dis[j])
			 {
				 temp1= O_WAT_IN_C_dis[i];
				  O_WAT_IN_C_dis[i]= O_WAT_IN_C_dis[j];
				  O_WAT_IN_C_dis[j]=temp1;
				  
				 
				 temp2 = O_WAT_IN_C_id[i];
				 O_WAT_IN_C_id[i] = O_WAT_IN_C_id[j];
				 O_WAT_IN_C_id[j]=temp2;
				 
				 temp3 = O_WAT_IN_C_coor[i];
				 O_WAT_IN_C_coor[i] = O_WAT_IN_C_coor[j];
				 O_WAT_IN_C_coor[j] = temp3;
				 
				 temp3 = H1_WAT_IN_C_coor[i];
				 H1_WAT_IN_C_coor[i]=H1_WAT_IN_C_coor[j];
				 H1_WAT_IN_C_coor[j]=temp3;
				 
				 temp3 = H2_WAT_IN_C_coor[i];
				 H2_WAT_IN_C_coor[i]=H2_WAT_IN_C_coor[j];
				 H2_WAT_IN_C_coor[j]=temp3;
			 }
		 }
	 }
	
	
	std::ofstream outfile1,outfile2;
	outfile1.open("wat_pdb");

    //~ std::cout << "number: " << O_H1_H2_WAT_IN_C_id[0].size() << std::endl;
	//~ for(index i = 0; i != O_H1_H2_WAT_IN_C_id[0].size(); ++i)
	for(index i = 0; i != num_wat; ++i)	
	{
		outfile1 <<"ATOM"<<std::setw(7) <<3*i+1<<std::setw(5)<<"O"<<std::setw(7)<<"WAT"<<std::setw(7)<<i+1<<std::setw(15)
		             <<O_WAT_IN_C_coor[i][0] << std::setw(15)<<O_WAT_IN_C_coor[i][1]<<std::setw(15)<<O_WAT_IN_C_coor[i][2]<<std::endl;
		outfile1 <<"ATOM"<<std::setw(7) <<3*i+2<<std::setw(5)<<"H1"<<std::setw(7)<<"WAT"<<std::setw(7)<<i+1<<std::setw(15)
		             <<H1_WAT_IN_C_coor[i][0] << std::setw(15)<<H1_WAT_IN_C_coor[i][1]<<std::setw(15)<<H1_WAT_IN_C_coor[i][2]<<std::endl;
		outfile1 <<"ATOM"<<std::setw(7) <<3*i+3<<std::setw(5)<<"H2"<<std::setw(7)<<"WAT"<<std::setw(7)<<i+1<<std::setw(15)
		             <<H2_WAT_IN_C_coor[i][0] << std::setw(15)<<H2_WAT_IN_C_coor[i][1]<<std::setw(15)<<H2_WAT_IN_C_coor[i][2]<<std::endl;		             		             
		outfile1 << "TER" <<std::endl;
		//~ std::cout <<std::setw(15) <<O_H1_H2_WAT_IN_C_id[0][i]<<std::setw(15)<<O_H1_H2_WAT_IN_C_id[1][i]<<std::setw(15)<<O_H1_H2_WAT_IN_C_id[2][i]<< std::endl;
		//~ std::cout <<std::setw(15) <<O_WAT_IN_C_coor[i][0]<<std::setw(15)<<O_WAT_IN_C_coor[i][1]<<std::setw(15)<<O_WAT_IN_C_coor[i][2]<< std::endl;		
		//~ std::cout <<std::setw(15) <<O_WAT_IN_C_id[i]<< std::endl;		
     }
     outfile1.close(); 
     
     outfile2.open("wat_c_pdb");
     for(index i = 0; i != 999; ++i)
     {
		outfile2 <<"ATOM"<<std::setw(7) <<i+1<<std::setw(5)<<"A"<<i+1<<std::setw(7)<<"gra"<<std::setw(7)<<1<<std::setw(15)
		             <<nc_data.atom_coordinate(frame, graphene_id[i])[0] << std::setw(15)<<nc_data.atom_coordinate(frame, graphene_id[i])[1]<<
		             std::setw(15)<<nc_data.atom_coordinate(frame, graphene_id[i])[2]<<std::endl;
	 }
	 for(index i = 999; i != 1400; ++i)
     {
		outfile2 <<"ATOM"<<std::setw(7) <<i+1<<std::setw(5)<<"B"<<i-999<<std::setw(7)<<"gra"<<std::setw(7)<<1<<std::setw(15)
		             <<nc_data.atom_coordinate(frame, graphene_id[i])[0] << std::setw(15)<<nc_data.atom_coordinate(frame, graphene_id[i])[1]<<
		             std::setw(15)<<nc_data.atom_coordinate(frame, graphene_id[i])[2]<<std::endl;
	 }
	 outfile2<<"TER"<<std::endl;
     for(index i = 1400; i != 2399; ++i)
     {
		outfile2 <<"ATOM"<<std::setw(7) <<i+1<<std::setw(5)<<"A"<<i -1399<<std::setw(7)<<"gra"<<std::setw(7)<<2<<std::setw(15)
		             <<nc_data.atom_coordinate(frame, graphene_id[i])[0] << std::setw(15)<<nc_data.atom_coordinate(frame, graphene_id[i])[1]<<
		             std::setw(15)<<nc_data.atom_coordinate(frame, graphene_id[i])[2]<<std::endl;
	 }
	 for(index i = 2399; i != 2800; ++i)
     {
		outfile2 <<"ATOM"<<std::setw(7) <<i+1<<std::setw(5)<<"B"<<i -2399<<std::setw(7)<<"gra"<<std::setw(7)<<2<<std::setw(15)
		             <<nc_data.atom_coordinate(frame, graphene_id[i])[0] << std::setw(15)<<nc_data.atom_coordinate(frame, graphene_id[i])[1]<<
		             std::setw(15)<<nc_data.atom_coordinate(frame, graphene_id[i])[2]<<std::endl;
	 }
	 outfile2<<"TER"<<std::endl;	 
	 
	 //~ std::cout<< "size:  "<<O_WAT_id.size() <<O_WAT_IN_C_id.size()<<std::endl;

     int sum = 0;
	 for(index i = 0; i != O_WAT_id.size(); ++i)
	 {
		 int judge = 0;
		 for(index j =0; j!= num_wat; ++j)
		 {
			 if(O_WAT_id[i] == O_WAT_IN_C_id[j])
			 {

				 judge = 1;   
			 }
		 }
		 if(judge == 0)
		 {
			 ++sum;
			  outfile2 <<"ATOM"<<std::setw(7) <<3*(sum-1)+2801<<std::setw(5)<<"O"<<std::setw(7)<<"WAT"<<std::setw(7)<<(sum-1)+3<<std::setw(15)
		       <<nc_data.atom_coordinate(frame, O_WAT_id[i])[0] << std::setw(15)<<nc_data.atom_coordinate(frame, O_WAT_id[i])[1]
		       <<std::setw(15)<<nc_data.atom_coordinate(frame, O_WAT_id[i])[2]<<std::endl;
		      outfile2 <<"ATOM"<<std::setw(7) <<3*(sum-1)+2802<<std::setw(5)<<"H1"<<std::setw(7)<<"WAT"<<std::setw(7)<<(sum-1)+3<<std::setw(15)
		       <<nc_data.atom_coordinate(frame, O_WAT_id[i]+1)[0] << std::setw(15)<<nc_data.atom_coordinate(frame, O_WAT_id[i]+1)[1]
		       <<std::setw(15)<<nc_data.atom_coordinate(frame, O_WAT_id[i]+1)[2]<<std::endl;
		      outfile2 <<"ATOM"<<std::setw(7) <<3*(sum-1)+2803<<std::setw(5)<<"H2"<<std::setw(7)<<"WAT"<<std::setw(7)<<(sum-1)+3<<std::setw(15)
		       <<nc_data.atom_coordinate(frame, O_WAT_id[i]+2)[0] << std::setw(15)<<nc_data.atom_coordinate(frame, O_WAT_id[i]+2)[1]
		       <<std::setw(15)<<nc_data.atom_coordinate(frame, O_WAT_id[i]+2)[2]<<std::endl;	       
		       outfile2 << "TER" <<std::endl;
		 }
	 }
	  std::cout << sum<<std::endl;
     outfile2.close();     
	return 0;
}





