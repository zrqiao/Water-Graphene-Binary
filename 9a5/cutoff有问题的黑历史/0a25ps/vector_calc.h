#ifndef __VECTOR_CALCULATOR__
#define __VECTOR_CALCULATOR__

#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <algorithm>
#include <numeric>
#include <functional>

namespace vector_calc{

typedef std::vector<double>::size_type index;
typedef std::vector<double>::difference_type d_value;

const double PI=3.14159265358979323846264338327950288419716939937510582;

template<typename T>
inline T points_distance(const std::vector<T>& point1,
	const std::vector<T>& point2);
template<typename T>
inline std::vector<T> get_vector(const std::vector<T>& point1,
	const std::vector<T>& point2);
template<typename T>
inline std::vector<T> vector_add(const std::vector<T>& vec1,
	const std::vector<T>& vec2);
template<typename T>
inline std::vector<T> vector_minus(const std::vector<T>& vec1,
	const std::vector<T>& vec2);
template<typename T>
inline void vector_increase(std::vector<T>& vec1,
	const std::vector<T>& vec2);
template<typename T>
inline std::vector<T> vector_multiply(const std::vector<T>& vec1,
	const std::vector<T>& vec2);
template<typename T>
inline T dot_product(const std::vector<T>& vec1,
	const std::vector<T>& vec2);
template<typename T>
inline T vector_module(const std::vector<T>& vec);
template<typename T>
inline T vector_angle(const std::vector<T>& vec1,
	const std::vector<T>& vec2);
template<typename T,typename Tt>
inline std::vector<T> get_mass_center(const std::vector<std::vector<T> >& coordinates,
	const std::vector<Tt>& mass);
template <typename T,typename Tt>
inline std::vector<T> scaling_vector(const std::vector<T>& vec,Tt ratio);
template <typename T,typename Tt>
inline std::vector<T> vector_divide(const std::vector<T>& vec,Tt ratio);
template <typename T>
inline std::vector<T> average_coordinate(const std::vector<std::vector<T> >& coordinates);

inline double min(const std::vector<double>& input_vector);//added by dragoon
inline double max(const std::vector<double>& input_vector);//added by dragoon
double angle_by_3points(std::vector<double>& vertex_coordinate, std::vector<double>& side_p1_coord, 
												std::vector<double>& side_p2_coord);//added by dragoon
double PI_value(); //added by dragoon

template <typename T,typename Tt>
inline void scaling_itself(std::vector<T>& vec,Tt ratio);
inline std::vector<double> cross_product(const std::vector<double>& vec1,
	const std::vector<double>& vec2);
template <typename T>
inline std::vector<T> midpoint(const std::vector<T>& point1,
	const std::vector<T>& point2);
template <typename T>
inline std::vector<T> normalize_vector(const std::vector<T>& vec);
inline double circumradius(const std::vector<double>& point1,
	const std::vector<double>& point2,const std::vector<double>& point3);
inline std::vector<double> circumcenter(const std::vector<double>& point1,
	const std::vector<double>& point2,const std::vector<double>& point3);
inline double circumcircle(const std::vector<double>& point1,
	const std::vector<double>& point2,const std::vector<double>& point3,
	std::vector<double>& circumcenter);
inline std::vector<double> matrix_multiply_vector(const std::vector<std::vector<double> >& mat,
	const std::vector<double>& vec);
inline std::vector<double> rotate_x(const std::vector<double>& vec,double angle);
inline std::vector<double> rotate_y(const std::vector<double>& vec,double angle);
inline std::vector<double> rotate_z(const std::vector<double>& vec,double angle);
inline std::vector<double>::size_type interaction_index(
	std::vector<double>::size_type i,std::vector<double>::size_type j);
inline std::vector<double>::size_type inter_self_index(
	std::vector<double>::size_type i,std::vector<double>::size_type j);



// 将相互作用的二维坐标转化为一维坐标(No-self-interaction)
inline index interaction_index(index i,index j)
{
	if(i==j)
	{
		std::cerr<<"ERROR! The two interaction indexes cannot be equal!"<<std::endl;
		exit(-1);
	}
	std::vector<double>::size_type a,b;
	if(i<j)
	{
		a=i;
		b=j;
	}
	else
	{
		a=j;
		b=i;
	}
	return (b-1)*b/2+a;
}

// 将相互作用的二维坐标转化为一维坐标(include sefl-interaction)
inline index inter_self_index(index i,index j)
{
	std::vector<double>::size_type a,b;
	if(i<j)
	{
		a=i;
		b=j;
	}
	else
	{
		a=j;
		b=i;
	}
	return (b+1)*b/2+a;
}

template<class T>
struct multiply_constant
{
	const T constant;
	multiply_constant(T _constant):constant(_constant){}
	T operator()(const T& x)const{return x*constant;}
};

template<class T>
struct divide_constant
{
	const T constant;
	divide_constant(T _constant):constant(_constant){}
	T operator()(const T& x)const{return x/constant;}
};

// 两点间距
template<typename T>
inline T points_distance(const std::vector<T>& point1,
	const std::vector<T>& point2)
{
	std::vector<T> vec(get_vector(point1,point2));
	return sqrt(dot_product(vec,vec));
}

// 得到两点间的向量
template<typename T>
inline std::vector<T> get_vector(const std::vector<T>& point1,
	const std::vector<T>& point2)
{
	std::vector<T> vec;
	std::transform(point2.begin(),point2.end(),point1.begin(),
		std::back_inserter(vec),std::minus<T>());
	return vec;
}

// 两向量相加
template<typename T>
inline std::vector<T> vector_add(const std::vector<T>& vec1,
	const std::vector<T>& vec2)
{
	std::vector<T> vec_fin;
	std::transform(vec1.begin(),vec1.end(),vec2.begin(),
		std::back_inserter(vec_fin),std::plus<T>());
	return vec_fin;
}

// 两向量相减
template<typename T>
inline std::vector<T> vector_minus(const std::vector<T>& vec1,
	const std::vector<T>& vec2)
{
	std::vector<T> vec_fin;
	std::transform(vec1.begin(),vec1.end(),vec2.begin(),
		std::back_inserter(vec_fin),std::minus<T>());
	return vec_fin;
}

// 将向量vec2加到vec1
template<typename T>
inline void vector_increase(std::vector<T>& vec1,
	const std::vector<T>& vec2)
{
	std::transform(vec1.begin(),vec1.end(),vec2.begin(),vec1.begin(),
		std::plus<T>());
}

template<typename T>
inline std::vector<T> vector_multiply(const std::vector<T>& vec1,
	const std::vector<T>& vec2)
{
	std::vector<T> answer;
	std::transform(vec1.begin(),vec1.end(),vec2.begin(),
		std::back_inserter(answer),std::multiplies<T>());
	return answer;
}

// 向量点乘
template<typename T>
inline T dot_product(const std::vector<T>& vec1,const std::vector<T>& vec2)
{
	return std::inner_product(vec1.begin(),vec1.end(),vec2.begin(),
		static_cast<T>(0));
}

// 向量模
template<typename T>
inline T vector_module(const std::vector<T>& vec)
{
	return sqrt(std::inner_product(vec.begin(),vec.end(),vec.begin(),
		static_cast<T>(0)));
}

// 向量夹角
template<typename T>
inline T vector_angle(const std::vector<T>& vec1,
	const std::vector<T>& vec2)
{
	return dot_product(vec1,vec2)/
		(vector_module(vec1)*vector_module(vec2));
}

// 质心
template<typename T,typename Tt>
inline std::vector<T> get_mass_center(const std::vector<std::vector<T> >& coordinates,
	const std::vector<Tt>& mass)
{
	if(coordinates.size()!=mass.size())
	{
		std::cerr<<"ERROR ! The number of coordinates and mass mismatch!"<<std::endl;
		exit(-1);
	}
	T sum_mass=std::accumulate(mass.begin(),mass.end(),static_cast<T>(0));
	std::vector<T> sum_coor(inner_product(coordinates.begin(),
		coordinates.end(),mass.begin(),
		std::vector<T>(coordinates.front().size(),static_cast<T>(0)),
		vector_add<T>,scaling_vector<T,T>));
	return vector_divide(sum_coor,sum_mass);
}

// 向量平均值
template <typename T>
inline std::vector<T> average_coordinate(const std::vector<std::vector<T> >& coordinates)
{
	std::vector<T> sum(std::accumulate(coordinates.begin(),
		coordinates.end(),std::vector<T>(coordinates.front().size(),static_cast<T>(0)),
		vector_add<T>));
	return vector_divide(sum,coordinates.size());
}

//求最小值 by dragoon
inline double min(const std::vector<double>& input_vector)
{
	double min_value=*(input_vector.begin());
	for(std::vector<double>::const_iterator input_it=input_vector.begin();input_it!=input_vector.end();++input_it)
	{
		if( (*input_it)<min_value ) min_value=*input_it;
	}
	return min_value;
}

//求最大值 by dragoon
inline double max(const std::vector<double>& input_vector)
{
	double max_value=*(input_vector.begin());
	for(std::vector<double>::const_iterator input_it=input_vector.begin();input_it!=input_vector.end();++input_it)
	{
		if( (*input_it)>max_value ) max_value=*input_it;
	}
	return max_value;
}

//3个点A-B-C形成的角BAC的角度(in degree) 用于计算氢键判据等 by dragoon
double angle_by_3points(std::vector<double>& vertex_coordinate, std::vector<double>& side_p1_coord, std::vector<double>& side_p2_coord)
{
	std::vector<double> lateral1=vector_minus(vertex_coordinate, side_p1_coord),
									lateral2=vector_minus(vertex_coordinate, side_p2_coord);
	double cos_theta=dot_product(lateral1,lateral2)/sqrt(vector_module(lateral1)*vector_module(lateral2));
	return 180.0*acos(cos_theta)/PI;
}

//返回PI by dragoon
double PI_value()
{
	return PI;
}

// 显示向量
template <typename T>
void show_vector(const std::vector<T>& vec)
{
	std::cout<<"["<<vec.size()<<"](";
	for(std::vector<double>::size_type i=0;i!=vec.size();++i)
	{
		if(i!=0)
			std::cout<<",";
		std::cout<<vec[i];
	}
	std::cout<<")";
}

// 向量乘以常数
template <typename T,typename Tt>
inline std::vector<T> scaling_vector(const std::vector<T>& vec,Tt ratio)
{
	std::vector<T> scaling;
	std::transform(vec.begin(),vec.end(),std::back_inserter(scaling),
		multiply_constant<T>(static_cast<T>(ratio)));
	//~ std::transform(vec.begin(),vec.end(),std::back_inserter(scaling),
		//~ [&](const T& x)->T{return x*ratio;});
	return scaling;
}

template <typename T,typename Tt>
inline std::vector<T> vector_divide(const std::vector<T>& vec,Tt ratio)
{
	std::vector<T> scaling;
	std::transform(vec.begin(),vec.end(),std::back_inserter(scaling),
		divide_constant<T>(static_cast<T>(ratio)));
	//~ std::transform(vec.begin(),vec.end(),std::back_inserter(scaling),
		//~ [&](const T& x)->T{return x/ratio;});
	return scaling;
}

// 向量本身乘以常数
template <typename T,typename Tt>
inline void scaling_itself(std::vector<T>& vec,Tt ratio)
{
	std::transform(vec.begin(),vec.end(),vec.begin(),
		multiply_constant<T>(static_cast<T>(ratio)));
	//~ std::transform(vec.begin(),vec.end(),vec.begin(),
		//~ [&](const T& x)->T{return x*ratio;});
}

// 向量叉乘
inline std::vector<double> cross_product(const std::vector<double>& vec1,
	const std::vector<double>& vec2)
{
	std::vector<double> product(3,0.0);
	product[0]=vec1[1]*vec2[2]-vec1[2]*vec2[1];
	product[1]=-vec1[0]*vec2[2]+vec1[2]*vec2[0];
	product[2]=vec1[0]*vec2[1]-vec1[1]*vec2[0];
	return product;
}

// 中点
template <typename T>
inline std::vector<T> midpoint(const std::vector<T>& point1,
	const std::vector<T>& point2)
{
	std::vector<T> sum(vector_add(point1,point2));
	return vector_divide(sum,2);
}

// 向量归一化
template <typename T>
inline std::vector<T> normalize_vector(const std::vector<T>& vec)
{
	return vector_divide(vec,vector_module(vec));
}

// 外接圆半径
inline double circumradius(const std::vector<double>& point1,
	const std::vector<double>& point2,const std::vector<double>& point3)
{
	std::vector<double> vertical_vec(cross_product(
		get_vector(point2,point1),get_vector(point3,point2)));
	double vertical_length=vector_module(vertical_vec);
	return points_distance(point1,point2)*points_distance(point2,point3)*
		points_distance(point3,point1)/vertical_length/2.0;
}

// 外接圆圆心
inline std::vector<double> circumcenter(const std::vector<double>& point1,
	const std::vector<double>& point2,const std::vector<double>& point3)
{
	std::vector<double> vertical_vec(cross_product(
		get_vector(point2,point1),get_vector(point3,point2)));
	double vertical_length=vector_module(vertical_vec);
	double c1=pow(points_distance(point2,point3),2)*
		dot_product(get_vector(point2,point1),get_vector(point3,point1))/
		(vertical_length*vertical_length)/2.0;
	double c2=pow(points_distance(point1,point3),2)*
		dot_product(get_vector(point1,point2),get_vector(point3,point2))/
		(vertical_length*vertical_length)/2.0;
	double c3=pow(points_distance(point1,point2),2)*
		dot_product(get_vector(point1,point3),get_vector(point2,point3))/
		(vertical_length*vertical_length)/2.0;
	std::vector<double> center;
	center.push_back(c1*point1[0]+c2*point2[0]+c3*point3[0]);
	center.push_back(c1*point1[1]+c2*point2[1]+c3*point3[1]);
	center.push_back(c1*point1[2]+c2*point2[2]+c3*point3[2]);
	return center;
}

// 外界元圆心和半径
inline double circumcircle(const std::vector<double>& point1,
	const std::vector<double>& point2,const std::vector<double>& point3,
	std::vector<double>& circumcenter)
{
	std::vector<double> vertical_vec(cross_product(
		get_vector(point2,point1),get_vector(point3,point2)));
	double vertical_length=vector_module(vertical_vec);
	double c1=pow(points_distance(point2,point3),2)*
		dot_product(get_vector(point2,point1),get_vector(point3,point1))/
		(vertical_length*vertical_length)/2.0;
	double c2=pow(points_distance(point1,point3),2)*
		dot_product(get_vector(point1,point2),get_vector(point3,point2))/
		(vertical_length*vertical_length)/2.0;
	double c3=pow(points_distance(point1,point2),2)*
		dot_product(get_vector(point1,point3),get_vector(point2,point3))/
		(vertical_length*vertical_length)/2.0;
	circumcenter.resize(0);
	circumcenter.push_back(c1*point1[0]+c2*point2[0]+c3*point3[0]);
	circumcenter.push_back(c1*point1[1]+c2*point2[1]+c3*point3[1]);
	circumcenter.push_back(c1*point1[2]+c2*point2[2]+c3*point3[2]);
	
	return points_distance(point1,point2)*points_distance(point2,point3)*
		points_distance(point3,point1)/vertical_length/2.0;
}

//矩阵乘以向量
inline std::vector<double> matrix_multiply_vector(const std::vector<std::vector<double> >& mat,
	const std::vector<double>& vec)
{
	std::vector<double> product;
	for(std::vector<double>::size_type i=0;i!=mat.size();++i)
		product.push_back(dot_product(mat[i],vec));
	return product;
}

//沿X轴旋转
inline std::vector<double> rotate_x(const std::vector<double>& vec,double angle)
{
	angle*=PI/180.0;
	std::vector<std::vector<double> > rotate_matrix(3,std::vector<double>(3,0));
	rotate_matrix[0][0]=1;
	rotate_matrix[1][1]=cos(angle);
	rotate_matrix[1][2]=-sin(angle);
	rotate_matrix[2][1]=sin(angle);
	rotate_matrix[2][2]=cos(angle);
	return matrix_multiply_vector(rotate_matrix,vec);
}

//沿Y轴旋转
inline std::vector<double> rotate_y(const std::vector<double>& vec,double angle)
{
	angle*=PI/180.0;
	std::vector<std::vector<double> > rotate_matrix(3,std::vector<double>(3,0));
	rotate_matrix[1][1]=1;
	rotate_matrix[0][0]=cos(angle);
	rotate_matrix[0][2]=-sin(angle);
	rotate_matrix[2][0]=sin(angle);
	rotate_matrix[2][2]=cos(angle);
	return matrix_multiply_vector(rotate_matrix,vec);
}

//沿Z轴旋转
inline std::vector<double> rotate_z(const std::vector<double>& vec,double angle)
{
	angle*=PI/180.0;
	std::vector<std::vector<double> > rotate_matrix(3,std::vector<double>(3,0));
	rotate_matrix[2][2]=1;
	rotate_matrix[0][0]=cos(angle);
	rotate_matrix[0][1]=-sin(angle);
	rotate_matrix[1][0]=sin(angle);
	rotate_matrix[1][1]=cos(angle);
	return matrix_multiply_vector(rotate_matrix,vec);
}

}

#endif
