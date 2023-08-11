#include "../include/Interpolator3D.h"
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

double Interpolator3D::pos_of_grid_point (XYZ xyz, int index, DataGenerationConfig& config)
{
    std::string grid_spacing;
    int n(0);
    double min(0), max(0);

    if (xyz == X) { grid_spacing = config.x_grid_spacing; n = config.n_x; min = config.x_min; max = config.x_max; }
    else if (xyz == Y) { grid_spacing = config.y_grid_spacing; n = config.n_y; min = config.y_min; max = config.y_max; }
    else if (xyz == Z) { grid_spacing = config.z_grid_spacing; n = config.n_z; min = config.z_min; max = config.z_max; }

    // log only really makes sense for values >=0
    if (grid_spacing == "linear")
    {
        return min+(max-min)*(double(index))/(double(n-1));
    }
    else if (config.x_grid_spacing == "log")
    {
        if (min == 0)
        {
            min = 1.0e-3*(max-min)/(double(n-1));
            return min*exp(log(max/min)*(double(index))/double(n-1))-min;
        }
        else
        {
            return min*exp(log(max/min)*(double(index))/double(n-1));
        }
    }
    else
    {
        return 0;
        exit(0);
    }
}


void Interpolator3D::set_grid (DataGenerationConfig& config)
{
    n_x = config.n_x;
    n_y = config.n_y;
    n_z = config.n_z;
    
    // initializing 3d vector and resizing it to size n_x,n_y,n_z
    x_pos.resize(n_x);
    y_pos.resize(n_y);
    z_pos.resize(n_z);

    for (unsigned int i=0; i<n_x; i++)
    {
        x_pos[i] = pos_of_grid_point(X,i,config);
    }
    for (unsigned int j=0; j<n_y; j++)
    {
        y_pos[j] = pos_of_grid_point(Y,j,config);
    }
    for (unsigned int k=0; k<n_z; k++)
    {
        z_pos[k] = pos_of_grid_point(Z,k,config);
    }
}


void Interpolator3D::prepare_data_array()
{
    for (auto& inner_vector : data_array)
    {
        for (auto& innermost_vector : inner_vector)
        {
            innermost_vector.clear();
        }
        inner_vector.clear();
    }
    data_array.clear();

    data_array.resize(n_x);
    for (unsigned int i=0; i<n_x; i++)
    {
       data_array[i].resize(n_y);
    }
    for (unsigned int i=0; i<n_x; i++)
    {
        for (unsigned int j=0; j<n_y; j++)
        {
            data_array[i][j].resize(n_z);
        }
    }
}


void Interpolator3D::export_data (std::string filepath) 
{
    std::cout << "Exporting data_array to file..." << std::endl;

    std::ofstream out;
    out.open(filepath);
    if (!out.is_open()) { std::cerr << "Could not open given file. Aborting" << std::endl; exit(0); }

    // printing n_x, n_y, n_z into the first line
    out << "#N " << n_x << " " << n_y << " " << n_z << std::endl;

    for (unsigned int i=0; i<n_x; i++)
    {
        for (unsigned int j=0; j<n_y; j++)
        {
            for (unsigned int k=0; k<n_z; k++)
            {
                out << std::setprecision(10) << i << " " << j << " " << k << " " << x_pos[i] << " " << y_pos[j] << " " << z_pos[k] << " " << data_array[i][j][k] << std::endl;
            }
            out << std::endl;
        }
        out << std::endl;
    }

    out.close();

    std::cout << "Finished exporting data_array to file" << std::endl;
}


void Interpolator3D::import_data (std::string filepath)
{
    std::cout << "Importing data_array from file..." << std::endl;

    std::ifstream in;
    in.open(filepath);
    if (!in.is_open()) { std::cerr << "Could not open given file. Aborting" << std::endl; exit(0); }

    std::string line;
    std::vector<std::string> line_vec;

    // getting n_x, n_y, n_z from the first line
    std::getline(in, line);
    std::istringstream ss(line);
    std::string element;
    while (std::getline(ss, element, ' '))
    {
        line_vec.push_back(element);
    }
    if (line_vec[0] != "#N") { std::cerr << "Wrong data_array format, use generate_data() to generate data_array in the desired format. Aborting!" << std::endl; exit(0); }

    int n_x = std::stoi(line_vec[1]);
    int n_y = std::stoi(line_vec[2]);
    int n_z = std::stoi(line_vec[3]);
    x_pos.resize(n_x);
    y_pos.resize(n_y);
    z_pos.resize(n_z);

    prepare_data_array();

    int i,j,k;
    while(std::getline(in, line))
    {
        if (line != "")
        {
            std::istringstream ss(line);
            std::string element;
            line_vec.clear();
            while (std::getline(ss, element, ' '))
            {
                line_vec.push_back(element);
            }

            i = std::stoi(line_vec[0]);
            j = std::stoi(line_vec[1]);
            k = std::stoi(line_vec[2]);

            x_pos[i] = std::stod(line_vec[3]); // these assignments are done unnecessarily often but I cant be bothered to do this smarter
            y_pos[j] = std::stod(line_vec[4]);
            z_pos[k] = std::stod(line_vec[5]);

            data_array[i][j][k] = std::stod(line_vec[6]);
        }
    }
    in.close();

    std::cout << "Finished importing data_array from file" << std::endl;
}


void Interpolator3D::generate_data (double func(double x, double y, double z), DataGenerationConfig& config, bool progress_monitor)
{    
    set_grid(config);
    prepare_data_array();

    // filling data_array with function values
    #pragma omp parallel for
    for (unsigned int i=0; i<n_x; i++)
    {
        for (unsigned int j=0; j<n_y; j++)
        {
            for (unsigned int k=0; k<n_z; k++)
            {
                data_array[i][j][k] = func(x_pos[i],y_pos[j],z_pos[k]);
            }
        }
        if (progress_monitor && (i % 10 == 0))
        {
            std::cout << i+1 << std::endl;
        }
    }
}


IndicesVec Interpolator3D::find_indices_of_closest_lower_data_point (double x, double y, double z)
{
    int i_0(0), j_0(0), k_0(0);
    int i_sign(1), j_sign(1), k_sign(1);
    
    if (x<x_pos[0])
    {
        i_sign = -1;
    }
    if (y<y_pos[0])
    {
        j_sign = -1;
    }
    if (z<z_pos[0])
    {
        k_sign = -1;
    }
    
    while (i_sign*x>i_sign*safe_get_x_pos(i_0+bool(i_sign+1)))
    {
        i_0 += i_sign;
    }
    while (j_sign*y>j_sign*safe_get_y_pos(j_0+bool(j_sign+1)))
    {
        j_0 += j_sign;
    }
    while (k_sign*z>k_sign*safe_get_z_pos(k_0+bool(k_sign+1)))
    {
        k_0 += k_sign;
    }
    
    IndicesVec indices_vec;

    indices_vec.i = i_0;
    indices_vec.j = j_0;
    indices_vec.k = k_0;

    return indices_vec;
}


double Interpolator3D::safe_get_x_pos (int i)
{
    if (i<0) return x_pos[0]+double(i);
    else if (i>int(n_x-1)) return x_pos[n_x-(unsigned int)1]+double(i-int(n_x+1));
    else return x_pos[(unsigned int)i];
}


double Interpolator3D::safe_get_y_pos (int i)
{
    if (i<0) return y_pos[0]+double(i)*1.0;
    else if (i>int(n_y-1)) return y_pos[n_y-(unsigned int)1]+double(i-int(n_y+1));
    else return y_pos[(unsigned int)i];
}


double Interpolator3D::safe_get_z_pos (int i)
{
    if (i<0) return z_pos[0]+double(i)*1.0;
    else if (i>int(n_z-1)) return z_pos[n_z-(unsigned int)1]+double(i-int(n_z+1));
    else return z_pos[(unsigned int)i];
}


double Interpolator3D::safe_get_data_point (int i, int j, int k)
{
    if (i<0) { i = 0; }
    else if (i>int(n_x-1)) { i = n_x-1; }

    if (j<0) { j = 0; }
    else if (j>int(n_y-1)) { j = n_y-1; }

    if (k<0) { k = 0; }
    else if (k>int(n_z-1)) { k = n_z-1; }

    return data_array[(unsigned int)i][(unsigned int)j][(unsigned int)k];
}


double Interpolator3D::unicubic_interpolate (double p[4], double t[4], double z)
{
    return p[1] + z*(p[2]-p[0])/(t[2]-t[0]) + z*z*(-3.0*p[1]+3.0*p[2]-2.0*(p[2]-p[0])/(t[2]-t[0])-(p[3]-p[1])/(t[3]-t[1])) + z*z*z*(2.0*p[1]-2.0*p[2]+(p[2]-p[0])/(t[2]-t[0])+(p[3]-p[1])/(t[3]-t[1]));
}


double Interpolator3D::bicubic_interpolate (double p[4][4], double t_y[4], double t_z[4], double y, double z)
{
	double unicubic_result[4];

	unicubic_result[0] = unicubic_interpolate(p[0], t_z, z);
	unicubic_result[1] = unicubic_interpolate(p[1], t_z, z);
	unicubic_result[2] = unicubic_interpolate(p[2], t_z, z);
	unicubic_result[3] = unicubic_interpolate(p[3], t_z, z);

	return unicubic_interpolate(unicubic_result, t_y, y);
}


double Interpolator3D::tricubic_interpolate (double p[4][4][4], double t_x[4], double t_y[4], double t_z[4], double x, double y, double z)
{
	double bicubic_result[4];

	bicubic_result[0] = bicubic_interpolate(p[0], t_y, t_z, y, z);
	bicubic_result[1] = bicubic_interpolate(p[1], t_y, t_z, y, z);
	bicubic_result[2] = bicubic_interpolate(p[2], t_y, t_z, y, z);
	bicubic_result[3] = bicubic_interpolate(p[3], t_y, t_z, y, z);

    return unicubic_interpolate(bicubic_result, t_x, x);
}


double Interpolator3D::tricubic_get_value (double x, double y, double z)
{
    IndicesVec indices_vec = find_indices_of_closest_lower_data_point(x,y,z);

    double p[4][4][4];
    for (int i=0; i<4; i++)
    {
        for (int j=0; j<4; j++)
        {
            for (int k=0; k<4; k++)
            {
                p[i][j][k] = safe_get_data_point(i+indices_vec.i-1,j+indices_vec.j-1,k+indices_vec.k-1);
            }
        }
    }

    double t_x[4];
    double t_y[4];
    double t_z[4];

    for (int i=0; i<4; i++)
    {
        t_x[i] = safe_get_x_pos(i+indices_vec.i-1);
        t_y[i] = safe_get_y_pos(i+indices_vec.j-1);
        t_z[i] = safe_get_z_pos(i+indices_vec.k-1);
    }

    x = (x-t_x[1])/(t_x[2]-t_x[1]);
    y = (y-t_y[1])/(t_y[2]-t_y[1]);
    z = (z-t_z[1])/(t_z[2]-t_z[1]);

    t_x[0] = (t_x[0]-t_x[1])/(t_x[2]-t_x[1]);
    t_x[3] = (t_x[3]-t_x[1])/(t_x[2]-t_x[1]);
    t_x[1] = 0.0;
    t_x[2] = 1.0;

    t_y[0] = (t_y[0]-t_y[1])/(t_y[2]-t_y[1]);
    t_y[3] = (t_y[3]-t_y[1])/(t_y[2]-t_y[1]);
    t_y[1] = 0.0;
    t_y[2] = 1.0;

    t_z[0] = (t_z[0]-t_z[1])/(t_z[2]-t_z[1]);
    t_z[3] = (t_z[3]-t_z[1])/(t_z[2]-t_z[1]);
    t_z[1] = 0.0;
    t_z[2] = 1.0;

    double tricubic_result = tricubic_interpolate(p,t_x,t_y,t_z,x,y,z);

    return tricubic_result;
}