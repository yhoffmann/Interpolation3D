#include "../include/Interpolator3D.h"
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>


double Interpolator3D::pos_of_grid_point (XYZ xyz, luint index, DataGenerationConfig& config)
{
    std::string grid_spacing;
    luint n(0);
    double min(0), max(0);

    if (xyz == X) { grid_spacing = config.x_grid_spacing; n = config.n_x; min = config.x_min; max = config.x_max; }
    else if (xyz == Y) { grid_spacing = config.y_grid_spacing; n = config.n_y; min = config.y_min; max = config.y_max; }
    else if (xyz == Z) { grid_spacing = config.z_grid_spacing; n = config.n_z; min = config.z_min; max = config.z_max; }

    // log only really makes sense for values > 0
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
    luint n_x = config.n_x;
    luint n_y = config.n_y;
    luint n_z = config.n_z;
    
    // initializing 3d vector and resizing it to size n_x,n_y,n_z
    x_pos.resize(n_x);
    y_pos.resize(n_y);
    z_pos.resize(n_z);

    for (luint i=0; i<n_x; i++)
    {
        x_pos[i] = pos_of_grid_point(X,i,config);
    }
    for (luint j=0; j<n_y; j++)
    {
        y_pos[j] = pos_of_grid_point(Y,j,config);
    }
    for (luint k=0; k<n_z; k++)
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

    data_array.resize(x_pos.size());
    for (luint i=0; i<x_pos.size(); i++)
    {
       data_array[i].resize(y_pos.size());
    }
    for (luint i=0; i<x_pos.size(); i++)
    {
        for (luint j=0; j<y_pos.size(); j++)
        {
            data_array[i][j].resize(z_pos.size());
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
    out << "#N " << x_pos.size() << " " << y_pos.size() << " " << z_pos.size() << std::endl;

    for (luint i=0; i<x_pos.size()-1; i++)
    {
        for (luint j=0; j<y_pos.size()-1; j++)
        {
            for (luint k=0; k<z_pos.size()-1; k++)
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

    luint n_x = std::stoi(line_vec[1]);
    luint n_y = std::stoi(line_vec[2]);
    luint n_z = std::stoi(line_vec[3]);
    x_pos.resize(n_x);
    y_pos.resize(n_y);
    z_pos.resize(n_z);

    prepare_data_array();

    luint i,j,k;
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

            line_vec.clear();
        }
    }
    in.close();

    std::cout << "Finished importing data_array from file" << std::endl;
}


void Interpolator3D::generate_data (double func(double x, double y, double z), DataGenerationConfig& config, bool progress_monitor)
{
    luint n_x = config.n_x;
    luint n_y = config.n_y;
    luint n_z = config.n_z;
    
    set_grid(config);
    prepare_data_array();

    // filling data_array with function values
    #pragma omp parallel for
    for (luint i=0; i<n_x; i++)
    {
        for (luint j=0; j<n_y; j++)
        {
            for (luint k=0; k<n_z; k++)
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


std::vector<int> Interpolator3D::find_indices_of_closest_smaller_data_point (double x, double y, double z)
{
    // TODO figure out how to deal with known data_array point hit or if this special case needs to be handled at all

    // finding corners of cuboid which x,y,z is inside of
    luint i_0(x_pos.size()-1), j_0(y_pos.size()-1), k_0(z_pos.size()-1);
    for (luint i=0; i<x_pos.size()-1; i++)
    {
        if (x_pos[i+1] > x)
        {
            i_0 = i;
            break;
        }
    }
    for (luint j=0; j<y_pos.size()-1; j++)
    {
        if (y_pos[j+1] > y)
        {
            j_0 = j;
            break;
        }
    }
    for (luint k=0; k<z_pos.size()-1; k++)
    {
        if (z_pos[k+1] > z)
        {
            k_0 = k;
            break;
        }
    }

    std::vector<int> ret;

    ret.push_back(i_0);
    ret.push_back(j_0);
    ret.push_back(k_0);

    return ret;
}


double Interpolator3D::get_clamped_pos (XYZ xyz, int i)
{
    if (xyz == X)
    {
        if (i<0) return x_pos[0];
        else if (i>x_pos[x_pos.size()-1]) return x_pos[x_pos.size()-1];
        else return x_pos[i];
    }
    if (xyz == Y)
    {
        if (i<0) return y_pos[0];
        else if (i>y_pos[y_pos.size()-1]) return y_pos[y_pos.size()-1];
        else return y_pos[i];
    }
    if (xyz == Z)
    {
        if (i<0) return z_pos[0];
        else if (i>z_pos[z_pos.size()-1]) return z_pos[z_pos.size()-1];
        else return z_pos[i];
    }
}


double Interpolator3D::safe_get_data_point (int i, int j, int k) // ints by design, these are intended to be <0 sometimes, may need to be made long int but realistically not
{
    if (i<0) { i = 0; }
    else if (i>x_pos.size()-1) { i = x_pos.size()-1; }

    if (j<0) { j = 0; }
    else if (j>y_pos.size()-1) { j = y_pos.size()-1; }

    if (k<0) { k = 0; }
    else if (k>z_pos.size()-1) { k = z_pos.size()-1; }

    return data_array[i][j][k];
}


double Interpolator3D::trilinear_get_value (double x, double y, double z)
{
    std::vector<int> indices_vec = find_indices_of_closest_smaller_data_point(x,y,z);

    luint i_0(indices_vec[0]), j_0(indices_vec[1]), k_0(indices_vec[2]);

    double x_0(x_pos[i_0]), x_1(x_pos[i_0+1]), y_0(y_pos[j_0]), y_1(y_pos[j_0+1]), z_0(z_pos[k_0]), z_1(z_pos[k_0+1]);

    // slopes
    double dx = (x-x_0)/(x_1-x_0);
    double dy = (y-y_0)/(y_1-y_0);
    double dz = (z-z_0)/(z_1-z_0);

    // interpolating along x
    double c_00 = safe_get_data_point(i_0,j_0,k_0) * (1.0-dx) + safe_get_data_point(i_0+1,j_0,k_0) * dx;
    double c_01 = safe_get_data_point(i_0,j_0,k_0+1) * (1.0-dx) + safe_get_data_point(i_0+1,j_0,k_0+1) * dx;
    double c_10 = safe_get_data_point(i_0,j_0+1,k_0) * (1.0-dx) + safe_get_data_point(i_0+1,j_0+1,k_0) * dx;
    double c_11 = safe_get_data_point(i_0,j_0+1,k_0+1) * (1.0-dx) + safe_get_data_point(i_0+1,j_0+1,k_0+1) * dx;

    // interpolating along y
    double c_0 = c_00 * (1.0-dy) + c_10 * dy;
    double c_1 = c_01 * (1.0-dy) + c_11 * dy;

    // interpolating along z, this is the interpolated value
    return c_0 * (1.0-dz) + c_1 * dz;
}


double Interpolator3D::unicubic_interpolate (double p[4], double z)
{
	return p[1] + 0.5 * z*(p[2] - p[0] + z*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + z*(3.0*(p[1] - p[2]) + p[3] - p[0])));
}


double Interpolator3D::bicubic_interpolate (double p[4][4], double y, double z)
{
	double unicubic_result[4];

	unicubic_result[0] = unicubic_interpolate(p[0], z);
	unicubic_result[1] = unicubic_interpolate(p[1], z);
	unicubic_result[2] = unicubic_interpolate(p[2], z);
	unicubic_result[3] = unicubic_interpolate(p[3], z);

	return unicubic_interpolate(unicubic_result, y);
}


double Interpolator3D::tricubic_interpolate (double p[4][4][4], double x, double y, double z)
{
	double bicubic_result[4];

	bicubic_result[0] = bicubic_interpolate(p[0], y, z);
	bicubic_result[1] = bicubic_interpolate(p[1], y, z);
	bicubic_result[2] = bicubic_interpolate(p[2], y, z);
	bicubic_result[3] = bicubic_interpolate(p[3], y, z);
    
	return unicubic_interpolate(bicubic_result, x);
}


double Interpolator3D::unicubic_interpolate_nonreg (double p[4], double t_z[4], double z)
{
    double t[4] = {t_z[0], t_z[1], t_z[2], t_z[3]};

    t[0] = (t[0]-t[1])/(t[2]-t[1]);
    t[3] = (t[3]-t[1])/(t[2]-t[1]);
    t[1] = 0.0;
    t[2] = 1.0;

	return p[1] + z*(p[2]-p[0])/(t[2]-t[0]) + z*z*(-3.0*p[1]+3.0*p[2]-2.0*(p[2]-p[0])/(t[2]-t[0])-(p[3]-p[1])/(t[3]-t[1])) + z*z*z*(2.0*p[1]-2.0*p[2]+(p[2]-p[0])/(t[2]-t[0])+(p[3]-p[1])/(t[3]-t[1]));
}


double Interpolator3D::bicubic_interpolate_nonreg (double p[4][4], double t_y[4], double t_z[4], double y, double z)
{
	double unicubic_result[4];

	unicubic_result[0] = unicubic_interpolate_nonreg(p[0], t_z, z);
	unicubic_result[1] = unicubic_interpolate_nonreg(p[1], t_z, z);
	unicubic_result[2] = unicubic_interpolate_nonreg(p[2], t_z, z);
	unicubic_result[3] = unicubic_interpolate_nonreg(p[3], t_z, z);

	return unicubic_interpolate_nonreg(unicubic_result, t_y, y);
}


double Interpolator3D::tricubic_interpolate_nonreg (double p[4][4][4], double t_x[4], double t_y[4], double t_z[4], double x, double y, double z)
{
	double bicubic_result[4];

	bicubic_result[0] = bicubic_interpolate_nonreg(p[0], t_y, t_z, y, z);
	bicubic_result[1] = bicubic_interpolate_nonreg(p[1], t_y, t_z, y, z);
	bicubic_result[2] = bicubic_interpolate_nonreg(p[2], t_y, t_z, y, z);
	bicubic_result[3] = bicubic_interpolate_nonreg(p[3], t_y, t_z, y, z);
    
	return unicubic_interpolate_nonreg(bicubic_result, t_x, x);
}


double Interpolator3D::bicubic_unilinear_get_value (double x, double y, double z)
{
    std::vector<int> indices_vec = find_indices_of_closest_smaller_data_point(x,y,z);

    luint i_0(indices_vec[0]), j_0(indices_vec[1]), k_0(indices_vec[2]);

    double x_0(x_pos[i_0]), x_1(x_pos[i_0+1]), y_0(y_pos[j_0]), y_1(y_pos[j_0+1]), z_0(z_pos[k_0]), z_1(z_pos[k_0+1]);

    double p_z_0[4][4];
    double p_z_1[4][4];
    for (int i=0; i<=3; i++)
    {
        for (int j=0; j<=3; j++)
        {
            p_z_0[i][j] = safe_get_data_point(i_0+i-1,j_0+j-1,k_0);
            p_z_1[i][j] = safe_get_data_point(i_0+i-1,j_0+j-1,k_0+1);
        }
    }

    x = (x-x_0)/(x_1-x_0);
    y = (y-y_0)/(y_1-y_0);
    z = (z-z_0)/(z_1-z_0);

    double bicubic_result_0 = bicubic_interpolate(p_z_0,x,y);
    double bicubic_result_1 = bicubic_interpolate(p_z_1,x,y);

    return bicubic_result_0+z*(bicubic_result_1-bicubic_result_0);
}


void Interpolator3D::setup_gsl_interp()
{
    gsl_3d.spline_1d = gsl_spline_alloc(gsl_3d.UNICUBIC,4);
    gsl_3d.spline_2d = gsl_spline2d_alloc(gsl_3d.BICUBIC,4,4);

    gsl_3d.x_acc = gsl_interp_accel_alloc();
    gsl_3d.y_acc = gsl_interp_accel_alloc();
    gsl_3d.z_acc = gsl_interp_accel_alloc();
}


void Interpolator3D::init_gsl_interp (int i_0, int j_0, int k_0)
{
    for (int i=i_0-1; i<i_0+3; i++)
    {
        gsl_3d.x_pos[i+1] = get_clamped_pos(X,i+i_0);
        gsl_3d.y_pos[i+1] = get_clamped_pos(Y,i+j_0);
    }

    for (int i=0; i<4; i++)
    {
        for (int j=0; j<4; j++)
        {
            gsl_3d.xy_array[i+4*j] = safe_get_data_point(i,j,k_0);
            std::cout << i << " " << j << " " << safe_get_data_point(i,j,k_0) << std::endl;
        }
    }

    gsl_spline2d_init(gsl_3d.spline_2d,gsl_3d.x_pos,gsl_3d.y_pos,gsl_3d.xy_array,4,4);
}


double Interpolator3D::tricubic_gsl_get_value (double x, double y, double z)
{
    std::vector<int> indices_vec = find_indices_of_closest_smaller_data_point(x,y,z);

    int i_0(indices_vec[0]), j_0(indices_vec[1]), k_0(indices_vec[2]);

    for (int k=0; k<4; k++)
    {
        init_gsl_interp(i_0,j_0,k_0);
        gsl_3d.z_array[k] = gsl_spline2d_eval(gsl_3d.spline_2d,x,y,gsl_3d.x_acc,gsl_3d.y_acc);
    }

    double ret = gsl_spline_eval(gsl_3d.spline_1d,z,gsl_3d.z_acc);

    return ret;
}


double Interpolator3D::tricubic_get_value (double x, double y, double z)
{
    std::vector<int> indices_vec = find_indices_of_closest_smaller_data_point(x,y,z);

    luint i_0(indices_vec[0]), j_0(indices_vec[1]), k_0(indices_vec[2]);

    double x_0(x_pos[i_0]), x_1(x_pos[i_0+1]), y_0(y_pos[j_0]), y_1(y_pos[j_0+1]), z_0(z_pos[k_0]), z_1(z_pos[k_0+1]);

    double p[4][4][4];
    for (int i=0; i<=3; i++)
    {
        for (int j=0; j<=3; j++)
        {
            for (int k=0; k<=3; k++)
            {
                p[i][j][k] = safe_get_data_point(i+i_0-1,j+j_0-1,k+k_0-1);
            }
        }
    }

    x = (x-x_0)/(x_1-x_0);
    y = (y-y_0)/(y_1-y_0);
    z = (z-z_0)/(z_1-z_0);

    double tricubic_result = tricubic_interpolate(p,x,y,z);

    return tricubic_result;
}


double Interpolator3D::tricubic_get_value_nonreg (double x, double y, double z)
{
    std::vector<int> indices_vec = find_indices_of_closest_smaller_data_point(x,y,z);

    int i_0(indices_vec[0]), j_0(indices_vec[1]), k_0(indices_vec[2]);

    double x_0(x_pos[i_0]), x_1(x_pos[i_0+1]), y_0(y_pos[j_0]), y_1(y_pos[j_0+1]), z_0(z_pos[k_0]), z_1(z_pos[k_0+1]);

    double p[4][4][4];
    for (int i=0; i<=3; i++)
    {
        for (int j=0; j<=3; j++)
        {
            for (int k=0; k<=3; k++)
            {
                p[i][j][k] = safe_get_data_point(i+i_0-1,j+j_0-1,k+k_0-1);
            }
        }
    }

    x = (x-x_0)/(x_1-x_0);
    y = (y-y_0)/(y_1-y_0);
    z = (z-z_0)/(z_1-z_0);

    double tricubic_result = tricubic_interpolate_nonreg(p,x,y,z);

    return tricubic_result;
}