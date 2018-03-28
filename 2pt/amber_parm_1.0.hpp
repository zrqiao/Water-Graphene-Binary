/**********************************************************************
 * Reading information of Amber topological files (.parm7 or .prmtop)
 * 
 *  Y. Isaac Yang
 *  College of Chemistry & Molecular Engineering
 *  Peking University
 *  Email: helloyesterday@pku.edu.cn
 * 
 * Zhuoran Dragoon Long Edition: several member functions added 
 * 02.07.2015 (mm.dd.yyyy)
 * longdragoonae@gmail.com
 * 
 **********************************************************************/

#ifndef __AMBER_PARAMETER__
#define __AMBER_PARAMETER__

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <set>
#include <iterator>
#include <map>

class amber_parm
{
public:
	typedef std::vector<double>::size_type index;
	typedef std::vector<double>::difference_type d_value;
	explicit amber_parm():if_empty(true),set_salt(false){}
	explicit amber_parm(const std::string& parm_file):name(parm_file),
		set_salt(false){open(name);}

	amber_parm& open(const std::string& parm_file);
	amber_parm& clear();
	amber_parm& close() {if_empty=true;return *this;}
	const std::string& file_name() const {return name;}
	bool empty() const {return if_empty;}

	index atoms_number() const {return natom;}
	index types_number() const {return ntypes;}
	index residues_number() const {return nres;}
	index if_box() const {return ifbox;}
	index solute_residues_number() const {return iptres;}
	index molecules_number() const {return nspm;}
	index solute_molecules_number() const {return nspsol;}
	index solute_atoms_number() const {return nsolute;}
	inline index index_by_molecule(index mol_id) const;

	const std::vector<index>& get_pointers() const {return pointers;}
	index get_pointers(index id) const {return pointers[id];}
	
	const std::vector<std::string>& atom_name() const {return igraph;}
	const std::string atom_name(index atom_id) const
		{return igraph[atom_id];}
	
	const std::vector<double>& atom_charge() const {return charge;}
	double atom_charge(index atom_id) const {return charge[atom_id];}
	
	const std::vector<double>& atom_mass() const {return amass;}
	double atom_mass(index atom_id) const {return amass[atom_id];}
	
	const std::vector<std::string>& atom_type() const {return isymbl;}
	std::string atom_type(index atom_id) const {return isymbl[atom_id];}
	
	const std::vector<std::string>& residue_label() const
		{return lbres;}
	std::string residue_label(index res_id) const {return lbres[res_id];}
	
	const std::vector<index>& residue_pointer() const {return ipres;}
	index residue_pointer(index res_id) const {return ipres[res_id];}
	
	const std::vector<index>& atoms_per_molecule() const {return nsp;}
	index atoms_per_molecule(index id) const {return nsp[id];}
	
	const std::vector<std::vector<index> >& atoms_by_residue() const
		{return res_atoms;}
	const std::vector<index>& atoms_by_residue(index res_id) const
		{return res_atoms[res_id];}

	std::vector<std::string> residue_atom(index res_index);
	const std::vector<index>& atom_resid() const {return ares;}
	index atom_resid(index id) const {return ares[id];}
	std::vector<double> residue_mass();
	std::vector<std::string> name_by_molecule(index mol_index)const;
	std::vector<double> mass_by_molecule(index mol_index)const;
	std::vector<double> molecule_mass();
	const std::vector<std::vector<double> >& LJ_ACOEF() const
		{return lj_acoef;}
	double LJ_ACOEF(index id1,index id2) const
		{return lj_acoef[id1][id2];}
		
	double LJ_BCOEF(index id1,index id2) const
		{return lj_bcoef[id1][id2];}
		
	const std::vector<std::string>& type_name() const {return iac_name;}
	std::string type_name(index type_id) const
		{return iac_name[type_id];}
	
	const std::vector<double>& types_mass() const {return iac_mass;}
	double types_mass(index type_id) const {return iac_mass[type_id];}
	
	//~ std::vector<double> types_charge() const {return iac_charge;}
	const std::vector<double>& types_radius() const {return iac_r0;}
	double types_radius(index type_id) const {return iac_r0[type_id];}
	
	const std::vector<double>& types_welldepth() const
		{return iac_epsilon;}
	double types_welldepth(index type_id) const
		{return iac_epsilon[type_id];}
	
	double type_mass(const std::string& type) const
		{return iac_mass[map_index(itype,type)];}
	//~ double type_charge(const std::string& type) const
		//~ {return iac_charge[map_index(itype,type)];}
	double type_radius(const std::string& type) const
		{return iac_r0[map_index(itype,type)];}
	double type_welldepth(const std::string& type) const
		{return iac_epsilon[map_index(itype,type)];}
		
	const std::vector<double>& atom_radius() const {return lj_r0;}
	double atom_radius(index atom_id) const {return lj_r0[atom_id];}
	
	const std::vector<double>& atom_welldepth() const
		{return lj_epsilon;}
	double atom_welldepth(index atom_id) const
		{return lj_epsilon[atom_id];}
	
	const std::map<std::string,index>& types_id(){return itype;}
	amber_parm& show_type(std::ostream& output)
		{do_show_type(output);return *this;}
	const amber_parm& show_type(std::ostream& output) const
		{do_show_type(output);return *this;}
	amber_parm& show_atom(std::ostream& output)
		{do_show_atom(output);return *this;}
	const amber_parm& show_atom(std::ostream& output) const
		{do_show_atom(output);return *this;}

	index id_by_name(const std::string& name) const;
	index id_by_name(const std::string& name,index res_id) const;
	
	
	std::vector<index> id_by_type(const std::string& type) const
		{return index_by_string(isymbl,type);}

	std::vector<index> id_by_residue(index res_id) const;
	std::vector<double> mass_by_residue(index res_id) const;
	std::vector<double> charge_by_residue(index res_id) const;
	
	std::vector<index> id_by_molecule(index mol_id);
	std::vector<std::vector<index> > ids_by_resname(
		std::vector<std::string>& resname) const;
	std::vector<index> ids_by_resname(const std::string& resname) const;
	std::vector<index> ids_by_resname_atomname(
		const std::string& resname, const std::string atomname) const;//added by dragoon
	std::vector<index> ids_by_resid_atomname(const std::string& atomname, index resid) const;//added by dragoon
	std::vector<index> ids_by_resid_atomname(
		const std::string& atomname, index resid_beg, index resid_end) const;//added by dragoon
	double mass_by_resname_atomname(const std::string& resname, const std::string& atomname) const; //added by dragoon

	amber_parm& estimate_salts();
	amber_parm& set_molecule(const std::vector<index>& nsp_out);
	index salts_number() const;
	
	const std::vector<std::string>& cation_atoms() const {return ication;}
	const std::string& cation_atoms(index cation_id) const {return ication[cation_id];}
	const std::vector<index>& cation_index() const {return cat_id;}
	const std::map<std::string,index>& cation_info() const {return tcation;}
	index cation_number() const {return ication.size();}
	index cation_number(const std::string& cation_type) const;
	
	const std::vector<std::string>& anion_atoms() const {return ianion;}
	const std::string& anion_atoms(index anion_id) const {return ianion[anion_id];}
	const std::vector<index>& anion_index() const {return ani_id;}
	const std::map<std::string,index>& anion_info() const {return tanion;}
	index anion_number() const {return ianion.size();}
	index anion_number(const std::string& anion_type) const;

private:
	std::string name;
	bool if_empty;

	// Origin data in amber parameter file:
	index natom;	// Total number of atoms
	index ntypes;	// Total number of atom types
	index nres;		// Total numberr of residues
	index ifbox;
	index iptres;	// Total number of solute residues
	index nspm;		// Total number of molecules
	index nspsol;	// Total number of solute molecules

	std::vector<index> pointers;
	std::vector<std::string> igraph;	// Atom name
	std::vector<double> charge;			// Atom charge
	std::vector<double> amass;			// Atom mass
	std::vector<index> iac;				// Atom type index
	std::vector<long> ico;				// Nonbonded parm index
	std::vector<std::string> lbres;		// Residue label
	std::vector<index> ipres;			// Residue pointers
	std::vector<double> cn1;			// Lennard Jones A coefficient
	std::vector<double> cn2;			// Lennard Jones B coefficient
	std::vector<std::string> isymbl;	// Amber atom type
	std::vector<index> solvate_pointers;
	std::vector<index> nsp;				// Atoms per molecule

	// Self_defined data
	bool set_salt;
	index nsolute;
	index nsalt;

	std::vector<std::string> iac_name;
	//~ std::vector<double> iac_charge;
	std::vector<double> iac_mass;
	std::vector<std::vector<double> > lj_acoef;
	std::vector<std::vector<double> > lj_bcoef;
	std::vector<double> lj_r0;
	std::vector<double> lj_epsilon;
	std::vector<double> iac_r0;
	std::vector<double> iac_epsilon;
	std::vector<index> ares;
	std::map<std::string,index> itype;	// ID of each atom type
	std::vector<std::vector<index> > res_atoms; // atoms in each residue;
	std::vector<std::string> ication;
	std::vector<std::string> ianion;
	std::vector<index> cat_id;
	std::vector<index> ani_id;
	std::map<std::string,index> tcation;
	std::map<std::string,index> tanion;

	// Built-in function
	inline std::string read_data(std::ifstream& ifile,
		std::vector<std::string>& array_string);
	inline std::string read_data(std::ifstream& ifile,
		std::vector<index>& array_index);
	inline std::string read_data(std::ifstream& ifile,
		std::vector<long>& array_long);
	inline std::string read_data(std::ifstream& ifile,
		std::vector<double>& array_double);
	amber_parm& calc_coef(const std::vector<double>& cn,
		std::vector<std::vector<double> >& coef);
	amber_parm& calc_vdw();
	amber_parm& get_type_data();
	amber_parm& set_solute(const std::vector<index>& sol_info);
	std::vector<index> index_by_string(
		const std::vector<std::string>& dou_data,
		const std::string& str_data) const;
	std::vector<double> double_by_residue(
		const std::vector<double>& dou_data,index res_id) const;
	amber_parm& set_residue();
	inline index map_index(const std::map<std::string,index>& map_type,
		const std::string& str) const;
	amber_parm& set_resatom();

	void do_show_type(std::ostream& output) const;
	void do_show_atom(std::ostream& output) const;
};

amber_parm& amber_parm::open(const std::string& parm_file)
{
	std::ifstream ifile(parm_file.c_str());
	if(!ifile)
	{
		std::cerr<<"ERROR in amber_parm.hpp! Can't open file: "<<
			parm_file<<std::endl;
		if_empty=true;
	}
	else
	{
		if_empty=false;
		std::string probe;
		ifile>>probe;
		while(!ifile.eof()&&probe!="%FLAG")
			ifile>>probe;
		ifile>>probe;
		while(!ifile.eof())
		{
			if(probe=="POINTERS")
				probe=read_data(ifile,pointers);
			else if(probe=="ATOM_NAME")
				probe=read_data(ifile,igraph);
			else if(probe=="CHARGE")
				probe=read_data(ifile,charge);
			else if(probe=="MASS")
				probe=read_data(ifile,amass);
			else if(probe=="ATOM_TYPE_INDEX")
				probe=read_data(ifile,iac);
			else if(probe=="NONBONDED_PARM_INDEX")
				probe=read_data(ifile,ico);
			else if(probe=="RESIDUE_LABEL")
				probe=read_data(ifile,lbres);
			else if(probe=="RESIDUE_POINTER")
				probe=read_data(ifile,ipres);
			else if(probe=="LENNARD_JONES_ACOEF")
				probe=read_data(ifile,cn1);
			else if(probe=="LENNARD_JONES_BCOEF")
				probe=read_data(ifile,cn2);
			else if(probe=="AMBER_ATOM_TYPE")
				probe=read_data(ifile,isymbl);
			else if(probe=="SOLVENT_POINTERS")
				probe=read_data(ifile,solvate_pointers);
			else if(probe=="ATOMS_PER_MOLECULE")
				probe=read_data(ifile,nsp);
			else
				ifile>>probe;
		}
		ifile.clear();
		ifile.close();
		natom=pointers[0];
		ntypes=pointers[1];
		if(igraph.size()!=natom||amass.size()!=natom||charge.size()!=natom||
			iac.size()!=natom)
		{
			std::cerr<<"ERROR in amber_parm.hpp! "<<
				"Number of atoms mismatch "<<
				"(number of natom/igraph/amass/iac):"<<std::endl<<natom<<
				" / "<<igraph.size()<<" / "<<amass.size()<<" / "<<
				iac.size()<<std::endl;
			exit(-1);
		}
		for(index i=0;i!=iac.size();++i)
			--iac[i];
		for(index i=0;i!=ipres.size();++i)
			--ipres[i];
		if(ipres.front()!=0)
		{
			std::cerr<<"ERROR in amber_parm.hpp! "<<
				"The first pointer of residue is not 1!"<<std::endl;
			exit(-1);
		}
		set_residue();
		calc_vdw();
		get_type_data();
		nres=pointers[11];
		if(lbres.size()!=nres||ipres.size()!=nres)
		{
			std::cerr<<"ERROR in amber_parm.hpp! "<<
				"Number of residues mismatch: (number of "<<
				"record/label/pointer)"<<std::endl<<
				nres<<" / "<<lbres.size()<<" / "<<ipres.size()<<std::endl;
			exit(-1);
		}
		ifbox=pointers[27];
		if(ifbox)
			set_solute(solvate_pointers);
	}
	set_resatom();
	return *this;
}


inline std::string amber_parm::read_data(std::ifstream& ifile,
	std::vector<std::string>& array_string)
{
	array_string.resize(0);
	std::string probe;
	while(probe!="%FORMAT(20a4)")	//读取无用的数据“%FORMAT(20a4)”
	{
		getline(ifile,probe);
		std::istringstream idata(probe);
		idata>>probe;
	}
	getline(ifile,probe);
	while(probe[0]!='%')
	{
		for(std::string::size_type i=0;i<probe.size();i+=4)
		{
			std::string word=probe.substr(i,4);
			std::istringstream delspace(word);	// 删除空格
			delspace>>word;
			array_string.push_back(word);
		}
		getline(ifile,probe);
	}
	std::istringstream idata(probe);
	idata>>probe;
	idata>>probe;
	return probe;
}

inline std::string amber_parm::read_data(std::ifstream& ifile,
	std::vector<index>& array_index)
{
	array_index.resize(0);
	std::string probe;
	ifile>>probe;	//读取无用的数据“%FORMAT(10I8)”
	ifile>>probe;
	while(probe!="%FLAG")
	{
		std::istringstream idata(probe);
		index data;
		idata>>data;
		array_index.push_back(data);
		ifile>>probe;
	}
	getline(ifile,probe);
	std::istringstream idata(probe);
	idata>>probe;
	return probe;
}

inline std::string amber_parm::read_data(std::ifstream& ifile,
	std::vector<long>& array_index)
{
	array_index.resize(0);
	std::string probe;
	ifile>>probe;	//读取无用的数据“%FORMAT(10I8)”
	ifile>>probe;
	while(probe!="%FLAG")
	{
		std::istringstream idata(probe);
		long data;
		idata>>data;
		array_index.push_back(data);
		ifile>>probe;
	}
	getline(ifile,probe);
	std::istringstream idata(probe);
	idata>>probe;
	return probe;
}

inline std::string amber_parm::read_data(std::ifstream& ifile,
	std::vector<double>& array_double)
{
	array_double.resize(0);
	std::string probe;
	ifile>>probe;	//读取无用的数据“%FORMAT(5E16.8)”
	ifile>>probe;
	while(probe!="%FLAG")
	{
		std::istringstream idata(probe);
		double data;
		idata>>data;
		array_double.push_back(data);
		ifile>>probe;
	}
	getline(ifile,probe);
	std::istringstream idata(probe);
	idata>>probe;
	return probe;
}

inline amber_parm::index amber_parm::map_index(
	const std::map<std::string,index>& map_type,
	const std::string& str) const
{
	std::map<std::string,index>::const_iterator it=map_type.find(str);
	if(it!=map_type.end())
		return it->second;
	else
	{
		std::cerr<<"ERROR in amber_parm.hpp: Can't find type: "<<str<<std::endl;
		return ntypes;
	}
}

amber_parm& amber_parm::calc_coef(const std::vector<double>& cn,
	std::vector<std::vector<double> >& coef)
{
	if(cn1.size()!=ntypes*(ntypes+1)/2)
	{
		std::cerr<<"ERROR in amber_parm.hpp! "<<
			"The number of Lennard Jones parameters mismatch: "<<
			ntypes*(ntypes+1)/2<<" / "<<cn1.size()<<std::endl;
		exit(-1);
	}
	coef.assign(ntypes,std::vector<double>(ntypes));
	index ia,ib,icn;
	for(index i=0;i!=coef.size();++i)
	{
		for(index j=0;j!=coef[i].size();++j)
		{
			if(i>j)
			{
				ia=i;
				ib=j;
			}
			else
			{
				ia=j;
				ib=i;
			}
			icn=ia*(ia+1)/2+ib;
			coef[i][j]=cn[icn];
		}
	}
	return *this;
}

amber_parm& amber_parm::calc_vdw()
{
	calc_coef(cn1,lj_acoef);
	calc_coef(cn2,lj_bcoef);
	iac_r0.resize(ntypes);
	iac_epsilon.resize(ntypes);
	for(index i=0;i!=iac_r0.size();++i)
	{
		if(fabs(lj_acoef[i][i])>1e-307)
		{
			iac_r0[i]=pow(2*lj_acoef[i][i]/lj_bcoef[i][i],1.0/6.0)/2.0;
			iac_epsilon[i]=lj_bcoef[i][i]/(2.0*pow(2*iac_r0[i],6));
		}
		else
		{
			iac_r0[i]=0;
			iac_epsilon[i]=0;
		}
	}
	lj_r0.resize(natom);
	lj_epsilon.resize(natom);
	for(index i=0;i!=natom;++i)
	{
		lj_r0[i]=iac_r0[iac[i]];
		lj_epsilon[i]=iac_epsilon[iac[i]];
	}
	return *this;
}

amber_parm& amber_parm::get_type_data()
{
	iac_name.resize(ntypes);
	iac_mass.resize(ntypes);
	//~ iac_charge.resize(ntypes);
	std::set<index> type_index;
	for(index i=0;i!=ntypes;++i)
		type_index.insert(i);
	for(index i=0;i!=iac.size()&&type_index.size()!=0;++i)
	{
		if(type_index.count(iac[i]))
		{
			iac_name[iac[i]]=isymbl[i];
			iac_mass[iac[i]]=amass[i];
			//~ iac_charge[iac[i]]=charge[i];
			type_index.erase(iac[i]);
		}
	}
	for(index i=0;i!=natom;++i)
		itype[isymbl[i]]=iac[i];
	return *this;
}

std::vector<std::string> amber_parm::residue_atom(index res_index)
{
	if(res_index>=nres)
	{
		std::cerr<<"ERROR in amber_parm.hpp! "<<
			"The residue index should less than"<<nres<<" !"<<std::endl;
		exit(-1);
	}
	index initial=ipres[res_index];
	index fin;
	
	if(res_index==nres-1)
		fin=natom;
	else
		fin=ipres[res_index+1];
	return std::vector<std::string>(igraph.begin()+initial,
		igraph.begin()+fin);
}

void amber_parm::do_show_type(std::ostream& output) const
{
	output<<std::left<<std::setw(8)<<"Index"<<std::setw(8)<<"Name"<<
		std::setw(10)<<"Mass"<<std::setw(10)<<"Charge"<<std::setw(10)<<
		"Radius"<<std::setw(10)<<"Welldepth"<<std::endl;
	for(index i=0;i!=ntypes;++i)
		output<<std::left<<std::setw(8)<<i<<std::setw(8)<<iac_name[i]<<
		std::setw(10)<<iac_mass[i]<<std::setw(10)<<iac_r0[i]<<
		std::setw(10)<<iac_epsilon[i]<<std::endl;
}

void amber_parm::do_show_atom(std::ostream& output) const
{
	output<<std::left<<std::setw(8)<<"Index"<<std::setw(8)<<"Name"<<
		std::setw(8)<<"Type"<<std::setw(10)<<"Mass"<<std::setw(10)<<
		"Charge"<<std::setw(10)<<"Radius"<<std::setw(10)<<"Welldepth"<<
		std::endl;
	for(index i=0;i!=natom;++i)
		output<<std::left<<std::setw(8)<<i<<std::setw(8)<<igraph[i]<<
		std::setw(8)<<isymbl[i]<<std::setw(10)<<amass[i]<<
		std::setw(10)<<charge[i]<<std::setw(10)<<lj_r0[i]<<
		std::setw(10)<<lj_epsilon[i]<<std::endl;
}

amber_parm& amber_parm::set_solute(const std::vector<index>& sol_info)
{
	if(sol_info.size()!=3)
	{
		std::cerr<<"ERROR in amber_parm.hpp! "<<
			"Solvate information break down!"<<std::endl;
		exit(-1);
	}
	iptres=sol_info[0];
	nspm=sol_info[1];
	nspsol=sol_info[2]-1;
	if(iptres>nres)
	{
		std::cerr<<"ERROR in amber_parm.hpp! "<<
			"Number of solute residues is too large!"<<std::endl;
		exit(-1);
	}
	if(iptres==nres)
		nsolute=natom;
	else
		nsolute=ipres[iptres];
	return *this;
}

std::vector<amber_parm::index> amber_parm::index_by_string(
	const std::vector<std::string>& class_data,
	const std::string& string_data) const
{
	std::vector<index> index_data;
	for(index i=0;i!=natom;++i)
	{
		if(class_data[i]==string_data)
			index_data.push_back(i);
	}
	return index_data;
}

std::vector<amber_parm::index> amber_parm::id_by_residue(index res_id)
	const
{
	if(res_id>=ipres.size())
	{
		std::cerr<<"ERROR in amber_parm.hpp! "<<
			"The residue ID is too large!"<<std::endl;
		exit(-1);
	}
	index id_initial=ipres[res_id];
	index id_final;
	if(res_id==ipres.size()-1)
		id_final=natom;
	else
		id_final=ipres[res_id+1];
	std::vector<index> atom_id;
	for(index i=id_initial;i!=id_final;++i)
		atom_id.push_back(i);
	return atom_id;
}

std::vector<double> amber_parm::double_by_residue(
	const std::vector<double>& double_data,index res_id) const
{
	index id_initial=ipres[res_id];
	index id_final;
	if(res_id==ipres.size()-1)
		id_final=natom;
	else
		id_final=ipres[res_id+1];
	std::vector<double> final_data;
	for(index i=id_initial;i!=id_final;++i)
		final_data.push_back(double_data[i]);
	return final_data;
}

std::vector<double> amber_parm::mass_by_residue(index res_id) const
{
	if(res_id>=ipres.size())
	{
		std::cerr<<"ERROR in amber_parm.hpp::"<<
			"mass_by_residue(index res_id):"<<
			std::endl<<"\tThe residue ID is too large!"<<std::endl;
		exit(-1);
	}
	return double_by_residue(amass,res_id);
}

amber_parm& amber_parm::estimate_salts()
{
	if(!ifbox)
	{
		std::cerr<<"ERROR in amber_parm.hpp::estimate_salts():"<<
			std::endl<<"\tIt have no box information!"<<std::endl;
		exit(-1);
	}
	nsalt=0;
	for(index i=0;i!=nsolute;++i)
	{
		std::string atype=isymbl[i];
		if(atype.size()>0)
		{
			char symbol=atype[atype.size()-1];
			if(symbol=='+')
			{
				ication.push_back(atype);
				cat_id.push_back(i);
				tcation[atype]=tcation[atype]+1;
				++nsalt;
			}
			if(symbol=='-')
			{
				ianion.push_back(atype);
				ani_id.push_back(i);
				tanion[atype]=tanion[atype]+1;
				++nsalt;
			}
		}
	}
	set_salt=true;
	return *this;
}

amber_parm::index amber_parm::salts_number() const
{
	if(!set_salt)
	{
		std::cerr<<std::endl<<"ERROR in amber_parm::salts_number():"<<
			std::endl<<"\tSalts information have not defined!"<<
			std::endl;
		exit(-1);
	}
	else
		return nsalt;
}

inline amber_parm::index amber_parm::index_by_molecule(index mol_index)
	const
{
	if(ifbox==0)
	{
		std::cerr<<"ERROR in amber_parm::index_by_molecule()!"<<
			std::endl<<"\tNo molecule information!"<<std::endl;
		exit(-1);
	}
	if(mol_index>=nspm)
	{
		std::cerr<<std::endl<<
			"ERROR in amber_parm::molecule_names(index mol_index):"<<
			std::endl<<"\tThe molecule index is larger than "<<
			"molecules number!"<<std::endl<<"\tnspm: "<<nspm<<
			"\tmol_index: "<<mol_index<<std::endl;
		exit(-1);
	}
	index initial=0;
	for(index i=0;i!=mol_index;++i)
		initial+=nsp[i];
	return initial;
}

std::vector<std::string> amber_parm::name_by_molecule(index mol_index)
	const
{
	if(ifbox==0)
	{
		std::cerr<<"ERROR in amber_parm::name_by_molecule()!"<<
			std::endl<<"\tNo molecule information!"<<std::endl;
		exit(-1);
	}
	index initial=index_by_molecule(mol_index);
	return std::vector<std::string>(igraph.begin()+initial,
		igraph.begin()+initial+nsp[mol_index]);
}

std::vector<double> amber_parm::mass_by_molecule(index mol_index) const
{
	if(ifbox==0)
	{
		std::cerr<<"ERROR in amber_parm::mass_by_molecule()!"<<
			std::endl<<"\tNo molecule information!"<<std::endl;
		exit(-1);
	}
	index initial=index_by_molecule(mol_index);
	return std::vector<double>(amass.begin()+initial,
		amass.begin()+initial+nsp[mol_index]);
}

amber_parm& amber_parm::set_molecule(const std::vector<index>& nsp_out)
{
	if(ifbox!=0)
	{
		std::cerr<<"ERROR in amber_parm.hpp: "<<
			"The parameter file include molucule infomation!"
			<<std::endl;
		exit(-1);
	}
	nspm=nsp_out.size();
	index atm_num=0;
	for(index i=0;i!=nsp_out.size();++i)
		atm_num+=nsp_out[i];
	if(atm_num!=natom)
	{
		std::cerr<<"ERROR in amber_parm.hpp:"<<
			"The number of atoms mismatch!"<<std::endl;
		exit(-1);
	}
	nsp=nsp_out;
	return *this;
}

std::vector<amber_parm::index> amber_parm::id_by_molecule(index mol_id)
{
	if(ifbox==0)
	{
		std::cerr<<"ERROR in amber_parm.hpp! No molecule information!"<<
			std::endl;
		exit(-1);
	}
	if(mol_id>=nsp.size())
	{
		std::cerr<<"ERROR in amber_parm.hpp: "<<
			"The moleclue id is larger than molecules number!"<<
			std::endl;
		exit(-1);
	}
	index begin_id=0;
	for(index i=0;i!=mol_id;++i)
		begin_id+=nsp[i];
	std::vector<index> idbm;
	for(index i=0;i!=nsp[mol_id];++i)
		idbm.push_back(begin_id+i);
	return idbm;
}

amber_parm& amber_parm::set_residue()
{
	index count=0;
	for(index i=0;i!=natom;++i)
	{
		if(i==ipres[count])
			++count;
		ares.push_back(count-1);
	}
	return *this;
}

amber_parm::index amber_parm::id_by_name(const std::string& name) const
{
	for(index i=0;i!=natom;++i)
	{
		if(igraph[i]==name)
			return i;
	}
	std::cerr<<"ERROR in amber_parm.hpp:"<<
		"Can't find atom named: "<<name<<std::endl;
	exit(-1);
}

amber_parm::index amber_parm::id_by_name(const std::string& name,
	index res_id) const
{
	if(res_id>=nres)
	{
		std::cerr<<"ERROR in amber_parm.hpp: "<<
			"The residue ID is larger than total residues!"<<std::endl;
		exit(-1);
	}
	for(index i=ipres[res_id];i!=(res_id==nres-1?natom:ipres[res_id+1]);
		++i)
	{
		if(igraph[i]==name)
			return i;
	}
	std::cerr<<"ERROR in amber_parm.hpp: Can't find atom named \""<<
		name<<"\" in residue \""<<res_id<<"\"."<<std::endl;
	exit(-1);
}

amber_parm& amber_parm::set_resatom()
{
	index res_id=0;
	std::vector<index>::const_iterator ires=ipres.begin();
	for(index i=0;i!=natom;++i)
	{
		if(ires!=ipres.end()&&i==*ires)
		{
			++ires;
			res_atoms.push_back(std::vector<index>());
			if(i!=0)
				++res_id;
		}
		res_atoms.back().push_back(i);
	}
	return *this;
}

std::vector<std::vector<amber_parm::index> > amber_parm::ids_by_resname(
	std::vector<std::string>& resname) const
{
	std::vector<std::vector<index> > res_ids;
	std::map<std::string,index> record;
	index id=0;
	for(index i=0;i!=res_atoms.size();++i)
	{
		index id_now;
		if(record.find(lbres[i])==record.end())
		{
			id_now=id;
			resname.push_back(lbres[i]);
			record[lbres[i]]=id++;
			res_ids.push_back(std::vector<index>());
		}
		else
			id_now=record[lbres[i]];
		std::copy(res_atoms[i].begin(),res_atoms[i].end(),
			std::back_inserter(res_ids[id_now]));
	}
	return res_ids;
}

std::vector<amber_parm::index> amber_parm::ids_by_resname(
	const std::string& resname) const
{
	std::vector<index> ids;
	for(index i=0;i!=res_atoms.size();++i)
	{
		if(lbres[i]==resname)
			std::copy(res_atoms[i].begin(),res_atoms[i].end(),
				std::back_inserter(ids));
	}
	return ids;
}

std::vector<amber_parm::index> amber_parm::ids_by_resname_atomname(
	const std::string& resname, const std::string atomname) const//by dragoon
{
	std::vector<index> ids;
	index status_label=0;
	for(index i=0;i!=res_atoms.size();++i)
	{
		if(lbres[i]==resname)
		{
			for(std::vector<index>::const_iterator a_in_res_it=res_atoms[i].begin(); a_in_res_it!=res_atoms[i].end();++a_in_res_it)
			{
				if(igraph[*a_in_res_it]==atomname) 
				{
					ids.push_back(*a_in_res_it);
					status_label=1;
				}
			}
		}
	}
	if(status_label) return ids;
	else
	{
		std::cerr<<"ERROR! Residue name "<<resname<<" or atom name "<<atomname<<" do not exist."<<std::endl;
		exit(-1);
	}
}

std::vector<amber_parm::index> amber_parm::ids_by_resid_atomname(
	const std::string& atomname, index resid) const//by dragoon
{
	if(resid>=ipres.size())
	{
		std::cerr<<"ERROR in amber_parm.hpp! "<<
						"The residue ID is too large!"<<std::endl;
		exit(-1);
	}
	index id_initial=ipres[resid];
	index id_final;
	if(resid==ipres.size()-1) id_final=natom;
	else id_final=ipres[resid+1];
	std::vector<index> atom_id;
	for(index i=id_initial;i!=id_final;++i) 
	{
		if(igraph[i]==atomname) atom_id.push_back(i);
	}
	return atom_id;
}

std::vector<amber_parm::index> amber_parm::ids_by_resid_atomname(
	const std::string& atomname, index resid_beg, index resid_end) const//by dragoon
{
	if(resid_beg>=ipres.size())
	{
		std::cerr<<"ERROR in amber_parm.hpp! "<<
						"The beginning residue ID is too large!"<<std::endl;
		exit(-1);
	}
	else if(resid_end>=ipres.size())
	{
		std::cerr<<"ERROR in amber_parm.hpp! "<<
						"The ending residue ID is too large!"<<std::endl;
		exit(-1);
	}
	std::vector<index> atom_id;
	for(index current_resid=resid_beg;current_resid!=resid_end+1;++current_resid)
	{
		index id_initial=ipres[current_resid];
		index id_final;
		if(current_resid==ipres.size()-1) id_final=natom;
		else id_final=ipres[current_resid+1];
		for(index i=id_initial;i!=id_final;++i) 
		{
			if(igraph[i]==atomname) atom_id.push_back(i);
		}
	}
	return atom_id;
}

double amber_parm::mass_by_resname_atomname(const std::string& resname, const std::string& atomname) const
{//dragoon 
	index status_label=0, first_id;
	for(index i=0;i!=res_atoms.size();++i)
	{
		if(lbres[i]==resname)
		{
			for(std::vector<index>::const_iterator a_in_res_it=res_atoms[i].begin(); a_in_res_it!=res_atoms[i].end();++a_in_res_it)
			{
				if(igraph[*a_in_res_it]==atomname) 
				{
					first_id=*a_in_res_it;
					status_label=1;
					break;
				}
			}
		}
	}
	if(status_label) return amass[first_id];
	else
	{
		std::cerr<<"ERROR! Residue name "<<resname<<" or atom name "<<atomname<<" do not exist."<<std::endl;
		exit(-1);
	}
}


amber_parm::index amber_parm::cation_number(const std::string& cation_type) const
{
	std::map<std::string,index>::const_iterator i=tcation.find(cation_type);
	if(i==tcation.end())
		return 0;
	else
		return i->second;
}

amber_parm::index amber_parm::anion_number(const std::string& anion_type) const
{
	std::map<std::string,index>::const_iterator i=tanion.find(anion_type);
	if(i==tanion.end())
		return 0;
	else
		return i->second;
}

#endif
