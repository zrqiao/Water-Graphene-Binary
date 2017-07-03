/**********************************************************************
 * Read information of Amber NetCDF trajactory file.				  *
 * Designed by Yi Isaac Yang										  *
 **********************************************************************/

#ifndef __AMBER_NETCDF__
#define __AMBER_NETCDF__

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <set>
#include <cmath>
#include <netcdfcpp.h>

class nctraj
{
public:
	typedef std::vector<double>::size_type index;
	typedef std::vector<double>::difference_type d_value;
	
	explicit nctraj();
	explicit nctraj(const std::string& name);
	const std::string& file_name() const {return fname;}

	index frames_number() const {return frame;}
	index atoms_number() const {return atom;}
	index spatials_number() const {return spatial;}
	float start_time() const {return start;}
	float time_interval() const {return dt;}

	float time_step(index step) const;
	std::vector<std::vector<double> > coordinates(index step) const;
	std::vector<double> atom_coordinate(index step,index atom_id) const;
	double spatial_coordinate(index step,index atom_id,
		index coor_id) const;
	double spatial_coordinate(index step,index atom_id,
		const std::string& coor_name) const;
	std::vector<std::vector<double> > coordinates_by_id(index step,
		const std::vector<index>& id_array) const;
	std::vector<std::vector<double> > coordinates_by_id(index step,
		index initial_index,index atom_nums) const;
		
	std::vector<std::vector<double> > velocities(index step) const;
	std::vector<double> atom_velocity(index step,index atom_id) const;
	double spatial_velocity(index step,index atom_id,
		index velo_id) const;
	double spatial_velocity(index step,index atom_id,
		const std::string& velo_name) const;
	std::vector<std::vector<double> > velocities_by_id(index step,
		const std::vector<index>& id_array) const;
	std::vector<std::vector<double> > velocities_by_id(index step,
		index initial_index,index atom_nums) const;

	std::vector<double> cell_lengths(index step) const;
	std::vector<double> cell_angles(index step) const;
	std::vector<std::vector<double> > period_box(index step)const;
	
	void show_dimensions();
	void show_attributes();
	
	void close(){file.close();}
	bool is_valid(){return file.is_valid();}
	void check_time() const;
	bool check_step(index step) const;
private:
	const std::string fname;
	NcFile file;

	std::set<std::string> dimensions;
	std::set<std::string> attributes;
	
	index frame;
	index spatial;
	index atom;
	index label;
	index cell_spatial;
	index cell_angular;

	const double zero_diff;
	float start;
	float dt;
	float velocities_scale_factor;

	std::string time_unit;
	std::string coordinates_unit;
	std::string velocities_units;
	std::string cell_lengths_unit;
	std::string cell_angles_unit;

	std::string title;
	std::string application;
	std::string program;
	std::string programVersion;
	std::string Conventions;
	std::string ConventionVersion;
	
	void get_dimension(index& dim,const std::string& dim_name);
	void get_attribute(std::string& att,const std::string& att_name);
};
nctraj::nctraj(const std::string& name):fname(name),
	file(name.c_str(),NcFile::ReadOnly),zero_diff(1e-7)
{
	if(!file.is_valid())
	{
		std::cerr<<"ERROR! Can't open file: "<<name<<std::endl;
		exit(-1);
	}
	
	for(index i=0;i!=static_cast<index>(file.num_dims());++i)
		dimensions.insert(file.get_dim(i)->name());
	for(index i=0;i!=static_cast<index>(file.num_atts());++i)
		attributes.insert(file.get_att(i)->name());

	get_dimension(frame,"frame");
	get_dimension(spatial,"spatial");
	get_dimension(atom,"atom");
	get_dimension(label,"label");
	get_dimension(cell_spatial,"cell_spatial");
	get_dimension(cell_angular,"cell_angular");

	get_attribute(title,"title");
	get_attribute(application,"application");
	get_attribute(program,"program");
	get_attribute(programVersion,"programVersion");
	get_attribute(Conventions,"Conventions");
	get_attribute(ConventionVersion,"ConventionVersion");
	
	const long int zero=0;
	float time12[2];
	file.get_var("time")->set_cur(zero);
	file.get_var("time")->get(time12,2);
	start=time12[0];
	dt=time12[1]-time12[0];
}

void nctraj::get_dimension(index& dim,const std::string& dim_name)
{
	if(dimensions.count(dim_name))
		dim=file.get_dim(dim_name.c_str())->size();
	else
		dim=0;
}

void nctraj::get_attribute(std::string& att,const std::string& att_name)
{
	if(attributes.count(att_name))
		att=file.get_att(att_name.c_str())->as_string(0);
	else
		att="No Data";
}

void nctraj::show_dimensions()
{
	std::cout<<std::endl<<"Dimensions in file: "<<fname<<std::endl;
	std::cout<<"-------------------------------------------"<<std::endl;
	std::cout<<"frame:\t\t"<<frame<<std::endl;
	std::cout<<"spatial:\t"<<spatial<<std::endl;
	std::cout<<"atom:\t\t"<<atom<<std::endl;
	std::cout<<"label:\t\t"<<label<<std::endl;
	if(cell_spatial==0||cell_angular==0)
	{
		std::cout<<"cell_spatial:\tNo Cell Information"<<std::endl;
		std::cout<<"cell_angular:\tNo Cell Information"<<std::endl;
	}
	else
	{
		std::cout<<"cell_spatial:\t"<<cell_spatial<<std::endl;
		std::cout<<"cell_angular:\t"<<cell_angular<<std::endl;
	}
	std::cout<<"-------------------------------------------"<<std::endl;
}

void nctraj::show_attributes()
{
	std::cout<<std::endl<<"Attributes in file: "<<fname<<std::endl;
	std::cout<<"-------------------------------------------"<<std::endl;
	std::cout<<"title:\t\t\t"<<title<<std::endl;
	std::cout<<"application:\t\t"<<application<<std::endl;
	std::cout<<"program:\t\t"<<program<<std::endl;
	std::cout<<"programVersion:\t\t"<<programVersion<<std::endl;
	std::cout<<"Conventions:\t\t"<<Conventions<<std::endl;
	std::cout<<"ConventionVersion:\t"<<ConventionVersion<<std::endl;
	std::cout<<"-------------------------------------------"<<std::endl;
}

std::vector<std::vector<double> > nctraj::coordinates(index step) const
{
	if(step>=frame)
	{
		std::cerr<<"ERROR! The step is larger than the total frame: "<<
			frame<<" : "<<step<<std::endl;
		exit(-1);
	}
	std::vector<std::vector<double> > coor_vec;
	for(index i=0;i!=atom;++i)
	{
		double coor_array[spatial];
		file.get_var("coordinates")->set_cur(step,i,0);
		file.get_var("coordinates")->get(coor_array,1,1,spatial);
		coor_vec.push_back(std::vector<double>(coor_array,coor_array+spatial));
	}
	return coor_vec;
}

std::vector<double> nctraj::atom_coordinate(index step,index atom_id) const
{
	if(step>=frame)
	{
		std::cerr<<"ERROR! The step is larger than the total frame: "<<
			frame<<" : "<<step<<std::endl;
		exit(-1);
	}
	if(atom_id>=atom)
	{
		std::cerr<<"ERROR! The atom ID is larger than the number of atoms: "<<
			atom<<" : "<<atom_id<<std::endl;
		exit(-1);
	}
	double coor_array[spatial];
	file.get_var("coordinates")->set_cur(step,atom_id,0);
	file.get_var("coordinates")->get(coor_array,1,1,spatial);
	std::vector<double> coor_vec(coor_array,coor_array+spatial);
	return coor_vec;
}

double nctraj::spatial_coordinate(index step,index atom_id,index coor_id) const
{
	if(step>=frame)
	{
		std::cerr<<"ERROR! The step is larger than the total frame: "<<
			frame<<" : "<<step<<std::endl;
		exit(-1);
	}
	if(atom_id>=atom)
	{
		std::cerr<<"ERROR! The atom ID is larger than the number of atoms: "<<
			atom<<" : "<<atom_id<<std::endl;
		exit(-1);
	}
	if(coor_id>=spatial)
	{
		std::cerr<<"ERROR! The coordiante ID is larger than spatial: "<<
			spatial<<" : "<<coor_id<<std::endl;
		exit(-1);
	}
	double coor;
	file.get_var("coordinates")->set_cur(step,atom_id,coor_id);
	file.get_var("coordinates")->get(&coor,1,1,1);
	return coor;
}

double nctraj::spatial_coordinate(index step,index atom_id,
	const std::string& coor_name) const
{
	if(step>=frame)
	{
		std::cerr<<"ERROR! The step is larger than the total frame: "<<
			frame<<" : "<<step<<std::endl;
		exit(-1);
	}
	if(atom_id>=atom)
	{
		std::cerr<<"ERROR! The atom ID is larger than the number of atoms: "<<
			atom<<" : "<<atom_id<<std::endl;
		exit(-1);
	}
	index coor_id;
	if(coor_name=="x"||coor_name=="X")
		coor_id=0;
	else if(coor_name=="y"||coor_name=="Y")
		coor_id=1;
	else if(coor_name=="z"||coor_name=="Z")
		coor_id=2;
	else
	{
		std::cerr<<"ERROR! The coordiantes name must be \"xyz\" or \"XYZ\"!"<<std::endl;
		exit(-1);
	}
	double coor;
	file.get_var("coordinates")->set_cur(step,atom_id,coor_id);
	file.get_var("coordinates")->get(&coor,1,1,1);
	return coor;
}

std::vector<std::vector<double> > nctraj::coordinates_by_id(index step,
	const std::vector<index>& id_array) const
{
	if(step>=frame)
	{
		std::cerr<<"ERROR! The step is larger than the total frame: "<<
			frame<<" : "<<step<<std::endl;
		exit(-1);
	}
	std::vector<std::vector<double> > coor;
	for(index i=0;i!=id_array.size();++i)
	{
		if(id_array[i]>=atom)
		{
			std::cerr<<"ERROR! The atom ID is larger than the number of atoms: "<<
				atom<<" : "<<id_array[i]<<std::endl;
			exit(-1);
		}
		coor.push_back(atom_coordinate(step,id_array[i]));
	}
	return coor;
}

std::vector<std::vector<double> > nctraj::coordinates_by_id(index step,
	index initial_index,index atom_nums) const
{
	if(step>=frame)
	{
		std::cerr<<"ERROR! The step is larger than the total frame: "<<
			frame<<" : "<<step<<std::endl;
		exit(-1);
	}
	if(initial_index>=atom)
	{
		std::cerr<<"ERROR! The atom ID is larger than the number of atoms: "<<
				atom<<" : "<<initial_index<<std::endl;
		exit(-1);
	}
	if(initial_index+atom_nums>atom)
	{
		std::cerr<<"ERROR! The atoms number is too large! "<<
				atom<<" : "<<initial_index+atom_nums<<std::endl;
		exit(-1);
	}
	std::vector<std::vector<double> > coor;
	for(index i=initial_index;i!=initial_index+atom_nums;++i)
		coor.push_back(atom_coordinate(step,i));
	return coor;
}
std::vector<std::vector<double> > nctraj::velocities(index step) const
{
	if(step>=frame)
	{
		std::cerr<<"ERROR! The step is larger than the total frame: "<<
			frame<<" : "<<step<<std::endl;
		exit(-1);
	}
	std::vector<std::vector<double> > velo;
	for(index i=0;i!=atom;++i)
	{
		double velo_array[spatial];
		file.get_var("velocities")->set_cur(step,i,0);
		file.get_var("velocities")->get(velo_array,1,1,spatial);
		velo.push_back(std::vector<double>(velo_array,velo_array+spatial));
	}
	return velo;
}

std::vector<double> nctraj::atom_velocity(index step,index atom_id) const
{
	if(step>=frame)
	{
		std::cerr<<"ERROR! The step is larger than the total frame: "<<
			frame<<" : "<<step<<std::endl;
		exit(-1);
	}
	if(atom_id>=atom)
	{
		std::cerr<<"ERROR! The atom ID is larger than the number of atoms: "<<
			atom<<" : "<<atom_id<<std::endl;
		exit(-1);
	}
	double velo_array[spatial];
	file.get_var("velocities")->set_cur(step,atom_id,0);
	file.get_var("velocities")->get(velo_array,1,1,spatial);
	std::vector<double> velo_vec(velo_array,velo_array+spatial);
	return velo_vec;
}

double nctraj::spatial_velocity(index step,index atom_id,index velo_id) const
{
	if(step>=frame)
	{
		std::cerr<<"ERROR! The step is larger than the total frame: "<<
			frame<<" : "<<step<<std::endl;
		exit(-1);
	}
	if(atom_id>=atom)
	{
		std::cerr<<"ERROR! The atom ID is larger than the number of atoms: "<<
			atom<<" : "<<atom_id<<std::endl;
		exit(-1);
	}
	if(velo_id>=spatial)
	{
		std::cerr<<"ERROR! The velocity ID is larger than spatial: "<<
			spatial<<" : "<<velo_id<<std::endl;
		exit(-1);
	}
	double velocity;
	file.get_var("velocities")->set_cur(step,atom_id,velo_id);
	file.get_var("velocities")->get(&velocity,1,1,1);
	return velocity;
}

double nctraj::spatial_velocity(index step,index atom_id,
	const std::string& velo_name) const
{
	if(step>=frame)
	{
		std::cerr<<"ERROR! The step is larger than the total frame: "<<
			frame<<" : "<<step<<std::endl;
		exit(-1);
	}
	if(atom_id>=atom)
	{
		std::cerr<<"ERROR! The atom ID is larger than the number of atoms: "<<
			atom<<" : "<<atom_id<<std::endl;
		exit(-1);
	}
	index velo_id;
	if(velo_name=="x"||velo_name=="X")
		velo_id=0;
	else if(velo_name=="y"||velo_name=="Y")
		velo_id=1;
	else if(velo_name=="z"||velo_name=="Z")
		velo_id=2;
	else
	{
		std::cerr<<"ERROR! The coordiantes name must be \"xyz\" or \"XYZ\"!"<<std::endl;
		exit(-1);
	}
	double velocity;
	file.get_var("velocities")->set_cur(step,atom_id,velo_id);
	file.get_var("velocities")->get(&velocity,1,1,1);
	return velocity;
}

std::vector<std::vector<double> > nctraj::velocities_by_id(index step,
	const std::vector<index>& id_array) const
{
	if(step>=frame)
	{
		std::cerr<<"ERROR! The step is larger than the total frame: "<<
			frame<<" : "<<step<<std::endl;
		exit(-1);
	}
	std::vector<std::vector<double> > velo;
	for(index i=0;i!=id_array.size();++i)
	{
		if(id_array[i]>=atom)
		{
			std::cerr<<"ERROR! The atom ID is larger than the number of atoms: "<<
				atom<<" : "<<id_array[i]<<std::endl;
			exit(-1);
		}
		velo.push_back(atom_velocity(step,id_array[i]));
	}
	return velo;
}

std::vector<std::vector<double> > nctraj::velocities_by_id(index step,
	index initial_index,index atom_nums) const
{
	if(step>=frame)
	{
		std::cerr<<"ERROR! The step is larger than the total frame: "<<
			frame<<" : "<<step<<std::endl;
		exit(-1);
	}
	if(initial_index>=atom)
	{
		std::cerr<<"ERROR! The atom ID is larger than the number of atoms: "<<
				atom<<" : "<<initial_index<<std::endl;
		exit(-1);
	}
	if(initial_index+atom_nums>atom)
	{
		std::cerr<<"ERROR! The atoms number is too large! "<<
				atom<<" : "<<initial_index+atom_nums<<std::endl;
		exit(-1);
	}
	std::vector<std::vector<double> > velo;
	for(index i=initial_index;i!=initial_index+atom_nums;++i)
		velo.push_back(atom_velocity(step,i));
	return velo;
}

std::vector<double> nctraj::cell_lengths(index step) const
{
	if(cell_spatial==0||cell_angular==0)
	{
		std::cerr<<"ERROR! This file have no cell information!"<<std::endl;
		exit(-1);
	}
	if(step>=frame)
	{
		std::cerr<<"ERROR! The step is larger than the total frame: "<<
			frame<<" : "<<step<<std::endl;
		exit(-1);
	}
	double box[cell_spatial];
	file.get_var("cell_lengths")->set_cur(step);
	file.get_var("cell_lengths")->get(box,1,cell_spatial);
	std::vector<double> cell_box(box,box+cell_spatial);
	return cell_box;
}

std::vector<double> nctraj::cell_angles(index step) const
{
	if(cell_spatial==0||cell_angular==0)
	{
		std::cerr<<"ERROR! This file have no cell information!"<<std::endl;
		exit(-1);
	}
	if(step>=frame)
	{
		std::cerr<<"ERROR! The step is larger than the total frame: "<<
			frame<<" : "<<step<<std::endl;
		exit(-1);
	}
	double box[cell_angular];
	file.get_var("cell_angles")->set_cur(step);
	file.get_var("cell_angles")->get(box,1,cell_angular);
	std::vector<double> box_angles(box,box+cell_angular);
	return box_angles;
}

std::vector<std::vector<double> > nctraj::period_box(index step) const
{
	const std::vector<double> box(cell_lengths(step));
	std::vector<std::vector<double> > period;
	std::vector<int> num(cell_spatial);
	num[0]=-1;
	num[1]=0;
	num[2]=1;
	for(index i=0;i!=cell_spatial;++i)
	{
		for(index j=0;j!=cell_spatial;++j)
		{
			for(index k=0;k!=cell_spatial;++k)
			{
				std::vector<double> tmp(cell_spatial);
				tmp[0]=num[i]*box[0];
				tmp[1]=num[j]*box[1];
				tmp[2]=num[k]*box[2];
				period.push_back(tmp);
			}
		}
	}
	return period;
}

float nctraj::time_step(index step) const
{
	if(step>=frame)
	{
		std::cerr<<"ERROR! The step is larger than the total frame: "<<
			frame<<" : "<<step<<std::endl;
		exit(-1);
	}
	float time_now;
	file.get_var("time")->set_cur(step);
	file.get_var("time")->get(&time_now,1);
	return time_now;
}

void nctraj::check_time() const
{
	float time_array[frame];
	const long int zero=0;
	file.get_var("time")->set_cur(zero);
	file.get_var("time")->get(time_array,frame);
	index ierror=0;
	std::vector<index> error_steps;
	for(index i=0;i!=frame;++i)
	{
		if(fabs(start+dt*i-time_array[i])>zero_diff)
		{
			error_steps.push_back(i);
			++ierror;
		}
	}
	if(!ierror)
		std::cout<<"\t"<<fname<<":\tCheck OK!"<<std::endl;
	else
	{
		std::string err_name=fname+".error_log";
		std::ofstream oerr(err_name.c_str());
		for(index i=0;i!=ierror;++i)
			oerr<<error_steps[i]<<std::endl;
		oerr.clear();
		oerr.close();
		std::cerr<<"Warning! File "<<fname<<" have "<<ierror<<
			" doubtful error snapshot(s)!"<<std::endl;
		std::cerr<<"   See more: "<<err_name<<std::endl;
	}
}

bool nctraj::check_step(index step) const
{
	float time_now;
	file.get_var("time")->set_cur(step);
	file.get_var("time")->get(&time_now,1);
	if(fabs(start+dt*step-time_now)<zero_diff)
		return true;
	else
		return false;
}

#endif
