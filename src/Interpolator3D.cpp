#include "../include/Interpolator3D.h"
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>


double Interpolator3D::inputs_for_pos(XYZ xyz, luint index, DataGenerationConfig& config)
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


void Interpolator3D::set_grid_pos (DataGenerationConfig& config)
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
        x_pos[i] = inputs_for_pos(X,i,config);
    }
    for (luint j=0; j<n_y; j++)
    {
        y_pos[j] = inputs_for_pos(Y,j,config);
    }
    for (luint k=0; k<n_z; k++)
    {
        z_pos[k] = inputs_for_pos(Z,k,config);
    }
}


void Interpolator3D::export_data_to_file (std::string filepath) 
{
    std::cout << "Exporting data to file" << std::endl;

    std::ofstream out;
    out.open(filepath);
    if (!out.is_open()) { std::cerr << "Could not open given file. Aborting" << std::endl; exit(0); }

    // printing n_x, n_y, n_z into the first line
    out << "#N " << x_pos.size() << " " << y_pos.size() << " " << z_pos.size() << std::endl;

    for (luint i=0; i<data.size(); i++)
    {
        for (luint j=0; j<data[0].size(); j++)
        {
            for (luint k=0; k<data[0][0].size(); k++)
            {
                out << std::setprecision(10) << i << " " << j << " " << k << " " << x_pos[i] << " " << y_pos[j] << " " << z_pos[k] << " " << data[i][j][k] << std::endl;
            }
            out << std::endl;
        }
        out << std::endl;
    }

    out.close();

    std::cout << "Finished exporting data to file" << std::endl;
}


void Interpolator3D::load_data (std::string filepath)
{
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
    if (line_vec[0] != "#N") { std::cerr << "Wrong data format, use generate_data() to generate data in the desired format. Aborting!" << std::endl; exit(0); }
    luint n_x = std::stoi(line_vec[1]);
    luint n_y = std::stoi(line_vec[2]);
    luint n_z = std::stoi(line_vec[3]);

    line_vec.clear();

    x_pos.resize(n_x);
    y_pos.resize(n_y);
    z_pos.resize(n_z);

    // getting the rest of the data
    vec_3d new_data;
    new_data.resize(n_x);
    for (luint i=0; i<n_x; i++)
    {
        new_data[i].resize(n_y);
    }
    for (luint i=0; i<n_x; i++)
    {
        for (luint j=0; j<n_y; j++)
        {
            new_data[i][j].resize(n_z);
        }
    }

    luint i,j,k;
    while(std::getline(in, line))
    {
        if (line != "")
        {
            std::istringstream ss(line);
            std::string element;
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

            new_data[i][j][k] = std::stod(line_vec[6]);

            line_vec.clear();
        }
    }
    in.close();

    data = new_data;
}


void Interpolator3D::generate_data (double func(double x, double y, double z), DataGenerationConfig& config, bool progress_monitor)
{
    luint n_x = config.n_x;
    luint n_y = config.n_y;
    luint n_z = config.n_z;
    
    // initializing 3d vector and resizing it to size n_x,n_y,n_z
    vec_3d new_data;
    new_data.resize(n_x);
    for (luint i=0; i<n_x; i++)
    {
        new_data[i].resize(n_y);
    }
    for (luint i=0; i<n_x; i++)
    {
        for (luint j=0; j<n_y; j++)
        {
            new_data[i][j].resize(n_z);
        }
    }

    set_grid_pos(config);

    // filling new_data with data from function func
    for (luint i=0; i<n_x; i++)
    {
        for (luint j=0; j<n_y; j++)
        {
            for (luint k=0; k<n_z; k++)
            {
                new_data[i][j][k] = func(x_pos[i],y_pos[j],z_pos[k]);
            }
        }
        if (progress_monitor && (i % 10 == 0))
        {
            std::cout << i+1 << std::endl;
        }
    }

    data = new_data;
}


double Interpolator3D::trilinear_get_value (double x, double y, double z)
{
    // TODO figure out how to deal with known data point hit or if this special case needs to be handled at all

    // handling inputs outside of known range
    if (x < x_pos[0] || x > x_pos[x_pos.size()-1]) { std::cout << "Requested value outside of known range! Aborting" << std::endl; exit(0); }
    if (y < y_pos[0] || y > y_pos[y_pos.size()-1]) { std::cout << "Requested value outside of known range! Aborting" << std::endl; exit(0); }
    if (z < z_pos[0] || z > z_pos[z_pos.size()-1]) { std::cout << "Requested value outside of known range! Aborting" << std::endl; exit(0); }

    // finding corners of cuboid which x,y,z is inside of
    double x_0(0), x_1(0), y_0(0), y_1(0), z_0(0), z_1(0);
    luint i_0(0), j_0(0), k_0(0);
    for (luint i=0; i<x_pos.size()-1; i++)
    {
        if (x_pos[i+1] > x)
        {
            x_0 = x_pos[i];
            x_1 = x_pos[i+1];
            i_0 = i;
            break;
        }
    }
    for (luint j=0; j<y_pos.size()-1; j++)
    {
        if (y_pos[j+1] > y)
        {
            y_0 = y_pos[j];
            y_1 = y_pos[j+1];
            j_0 = j;
            break;
        }
    }
    for (luint k=0; k<z_pos.size()-1; k++)
    {
        if (z_pos[k+1] > z)
        {
            z_0 = z_pos[k];
            z_1 = z_pos[k+1];
            k_0 = k;
            break;
        }
    }

    // slopes
    double dx = (x-x_0)/(x_1-x_0);
    double dy = (y-y_0)/(y_1-y_0);
    double dz = (z-z_0)/(z_1-z_0);

    // interpolating along x
    double c_00 = data[i_0][j_0][k_0] * (1.0-dx) + data[i_0+1][j_0][k_0] * dx;
    double c_01 = data[i_0][j_0][k_0+1] * (1.0-dx) + data[i_0+1][j_0][k_0+1] * dx;
    double c_10 = data[i_0][j_0+1][k_0] * (1.0-dx) + data[i_0+1][j_0+1][k_0] * dx;
    double c_11 = data[i_0][j_0+1][k_0+1] * (1.0-dx) + data[i_0+1][j_0+1][k_0+1] * dx;

    // interpolating along y
    double c_0 = c_00 * (1.0-dy) + c_10 * dy;
    double c_1 = c_01 * (1.0-dy) + c_11 * dy;

    // interpolating along z, this is the interpolated value
    return c_0 * (1.0-dz) + c_1 * dz;
}


double Interpolator3D::unicubic_interpolate (double p[4], double x)
{
	return p[1] + 0.5 * x*(p[2] - p[0] + x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + x*(3.0*(p[1] - p[2]) + p[3] - p[0])));
}

double Interpolator3D::bicubic_interpolate (double p[4][4], double x, double y)
{
	double unicubic_result[4];

	unicubic_result[0] = unicubic_interpolate(p[0], y);
	unicubic_result[1] = unicubic_interpolate(p[1], y);
	unicubic_result[2] = unicubic_interpolate(p[2], y);
	unicubic_result[3] = unicubic_interpolate(p[3], y);

	return unicubic_interpolate(unicubic_result, x);
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


double Interpolator3D::slope_at_vertex (XYZ t, int i, int j, int k)
{
    double dp, dt;
    double p_m1, p_1;
    luint index = (t==X) ? i : j;
    
    if (index==0)
    {
        if (t==X)
        {
            p_1 = data[i+1][j][k];
            p_m1 = data[i][j][k];
            dt = x_pos[i+1]-x_pos[i];
            return (p_1-p_m1)/dt;
        }
        else
        {
            p_1 = data[i][j+1][k];
            p_m1 = data[i][j][k];
            dt = y_pos[j+1]-y_pos[j];
            return (p_1-p_m1)/dt;
        }
    }
    else if (index==(t==X ? x_pos.size()-1 : y_pos.size()-1)) // TODO get rid of switch statements
    {
        switch (t)
        {
        case X:
            p_1 = data[i][j][k];
            p_m1 = data[i-1][j][k];
            dt = x_pos[i]-x_pos[i-1];
            return (p_1-p_m1)/dt;
        break;
        
        case Y:
            p_1 = data[i][j][k];
            p_m1 = data[i][j-1][k];
            dt = y_pos[j]-y_pos[j-1];
            return (p_1-p_m1)/dt;
        break;
        }
    }
    else
    {
        switch (t)
        {
        case X:
            p_1 = data[i+1][j][k];
            p_m1 = data[i-1][j][k];
            dt = x_pos[i+1]-x_pos[i-1];
            return (p_1-p_m1)/dt;
        break;
        
        case Y:
            p_1 = data[i][j+1][k];
            p_m1 = data[i][j-1][k];
            dt = y_pos[j+1]-y_pos[j-1];
            return (p_1-p_m1)/dt;
        break;
        }
    }
    return 0;
}




double Interpolator3D::slope_at_vertex (XYZ t_1, XYZ t_2, int i, int j, int k)
{
    double dp_x, dp_y, dt_x, dt_y;
    double p_m1m1, p_11, p_m11, p_1m1;

    if (i==0 && j==0)
    {
        p_11 = data[i+1][j+1][k];
        p_1m1 = data[i+1][j][k];
        p_m11 = data[i][j+1][k];
        p_m1m1 = data[i][j][k];

        dt_x = x_pos[i+1]-x_pos[i];
        dt_y = y_pos[j+1]-y_pos[j];

        return (p_11-p_1m1-p_m11+p_m1m1)/dt_x/dt_y;
    }
    else if (i==0 && j==y_pos.size()-1)
    {
        p_11 = data[i+1][j][k];
        p_1m1 = data[i+1][j-1][k];
        p_m11 = data[i][j][k];
        p_m1m1 = data[i][j-1][k];

        dt_x = x_pos[i+1]-x_pos[i];
        dt_y = y_pos[j]-y_pos[j-1];

        return (p_11-p_1m1-p_m11+p_m1m1)/dt_x/dt_y;
    }
    else if (i==x_pos.size()-1 && j==0)
    {
        p_11 = data[i][j+1][k];
        p_1m1 = data[i][j][k];
        p_m11 = data[i-1][j+1][k];
        p_m1m1 = data[i-1][j][k];

        dt_x = x_pos[i]-x_pos[i-1];
        dt_y = y_pos[j+1]-y_pos[j];

        return (p_11-p_1m1-p_m11+p_m1m1)/dt_x/dt_y;
    }
    else if (i==x_pos.size()-1 && j==y_pos.size()-1)
    {
        p_11 = data[i][j][k];
        p_1m1 = data[i][j-1][k];
        p_m11 = data[i-1][j][k];
        p_m1m1 = data[i-1][j-1][k];

        dt_x = x_pos[i]-x_pos[i-1];
        dt_y = y_pos[j]-y_pos[j-1];

        return (p_11-p_1m1-p_m11+p_m1m1)/dt_x/dt_y;
    }
    else if (i==0)
    {
        p_11 = data[i+1][j+1][k];
        p_1m1 = data[i+1][j-1][k];
        p_m11 = data[i][j+1][k];
        p_m1m1 = data[i][j-1][k];

        dt_x = x_pos[i+1]-x_pos[i];
        dt_y = y_pos[j+1]-y_pos[j-1];

        return (p_11-p_1m1-p_m11+p_m1m1)/dt_x/dt_y;
    }
    else if (j==0)
    {
        p_11 = data[i+1][j+1][k];
        p_1m1 = data[i+1][j][k];
        p_m11 = data[i-1][j+1][k];
        p_m1m1 = data[i-1][j][k];

        dt_x = x_pos[i+1]-x_pos[i-1];
        dt_y = y_pos[j+1]-y_pos[j];

        return (p_11-p_1m1-p_m11+p_m1m1)/dt_x/dt_y;
    }
    else
    {
        p_11 = data[i+1][j+1][k];
        p_1m1 = data[i+1][j-1][k];
        p_m11 = data[i-1][j+1][k];
        p_m1m1 = data[i-1][j-1][k];

        dt_x = x_pos[i+1]-x_pos[i-1];
        dt_y = y_pos[j+1]-y_pos[j-1];

        return (p_11-p_1m1-p_m11+p_m1m1)/dt_x/dt_y;
    }
}


double Interpolator3D::bicubic_get_value (double x, double y, luint k_0)
{
    // handling inputs outside of known range
    if (x < x_pos[0] || x > x_pos[x_pos.size()-1]) { std::cout << "Requested value outside of known range! Aborting" << std::endl; exit(0); }
    if (y < y_pos[0] || y > y_pos[y_pos.size()-1]) { std::cout << "Requested value outside of known range! Aborting" << std::endl; exit(0); }

    // finding corners of cuboid which x,y,z is inside of
    double x_0(0), x_1(0), y_0(0), y_1(0), z_0(0), z_1(0);
    luint i_0(0), j_0(0);
    for (luint i=0; i<x_pos.size()-1; i++)
    {
        if (x_pos[i+1] > x)
        {
            x_0 = x_pos[i];
            x_1 = x_pos[i+1];
            i_0 = i;
            break;
        }
    }
    for (luint j=0; j<y_pos.size()-1; j++)
    {
        if (y_pos[j+1] > y)
        {
            y_0 = y_pos[j];
            y_1 = y_pos[j+1];
            j_0 = j;
            break;
        }
    }

    double a_vec[16]{0.0};
    
    double dx = x_pos[i_0+1]-x_pos[i_0];
    double dy = y_pos[j_0+1]-y_pos[j_0];

    double x_vec[16] =
    {
        data[i_0][j_0][k_0],
        data[i_0+1][j_0][k_0],
        data[i_0][j_0+1][k_0],
        data[i_0+1][j_0+1][k_0],
        dx*slope_at_vertex(X,i_0,j_0,k_0),
        dx*slope_at_vertex(X,i_0+1,j_0,k_0),
        dx*slope_at_vertex(X,i_0,j_0+1,k_0),
        dx*slope_at_vertex(X,i_0+1,j_0+1,k_0),
        dy*slope_at_vertex(Y,i_0,j_0,k_0),
        dy*slope_at_vertex(Y,i_0+1,j_0,k_0),
        dy*slope_at_vertex(Y,i_0,j_0+1,k_0),
        dy*slope_at_vertex(Y,i_0+1,j_0+1,k_0),
        dx*dy*slope_at_vertex(X,Y,i_0,j_0,k_0),
        dx*dy*slope_at_vertex(X,Y,i_0+1,j_0,k_0),
        dx*dy*slope_at_vertex(X,Y,i_0,j_0+1,k_0),
        dx*dy*slope_at_vertex(X,Y,i_0+1,j_0+1,k_0)
    };

    for (int i=0; i<16; i++)
    {
        for (int j=0; j<16; j++)
        {
            a_vec[i] += BICUBIC_COEFFICIENTS[i][j] * x_vec[j];
        }
    }

    // normalizing the cubiod to a unit cube
    x = (x-x_0)/(x_1-x_0);
    y = (y-y_0)/(y_1-y_0);

    double interpolated_value(0);

    for (int i=0; i<=3; i++)
    {
        for (int j=0; j<=3; j++)
        {
            interpolated_value += a_vec[i+4*j] * std::pow(x,i) * std::pow(y,j);
        }
    }

    return interpolated_value;
}


double Interpolator3D::bicubic_unilinear_get_value (double x, double y, double z)
{
    if (z < z_pos[0] || z > z_pos[z_pos.size()-1]) { std::cout << "Requested value outside of known range! Aborting" << std::endl; exit(0); }
    double z_0(0), z_1(0);
    luint k_0(0);
    for (luint k=0; k<z_pos.size()-1; k++)
    {
        if (z_pos[k+1] > z)
        {
            z_0 = z_pos[k];
            z_1 = z_pos[k+1];
            k_0 = k;
            break;
        }
    }

    double bicubic_result_0 = bicubic_get_value(x,y,k_0);
    double bicubic_result_1 = bicubic_get_value(x,y,k_0+1);

    z = (z-z_0)/(z_pos[k_0+1]-z_pos[k_0]);

    return bicubic_result_0+z*(bicubic_result_1-bicubic_result_0);
}


double Interpolator3D::bicubic_unilinear_get_value_test (double x, double y, double z)
{
    // handling inputs outside of known range
    if (x < x_pos[0] || x > x_pos[x_pos.size()-1]) { std::cout << "Requested value outside of known range! Aborting" << std::endl; exit(0); }
    if (y < y_pos[0] || y > y_pos[y_pos.size()-1]) { std::cout << "Requested value outside of known range! Aborting" << std::endl; exit(0); }
    if (z < z_pos[0] || z > z_pos[z_pos.size()-1]) { std::cout << "Requested value outside of known range! Aborting" << std::endl; exit(0); }

    // finding corners of cuboid which x,y,z is inside of
    double x_0(0), x_1(0), y_0(0), y_1(0), z_0(0), z_1(0);
    luint i_0(0), j_0(0), k_0(0);
    for (luint i=0; i<x_pos.size()-1; i++)
    {
        if (x_pos[i+1] > x)
        {
            x_0 = x_pos[i];
            x_1 = x_pos[i+1];
            i_0 = i;
            break;
        }
    }
    for (luint j=0; j<y_pos.size()-1; j++)
    {
        if (y_pos[j+1] > y)
        {
            y_0 = y_pos[j];
            y_1 = y_pos[j+1];
            j_0 = j;
            break;
        }
    }
    for (luint k=0; k<z_pos.size()-1; k++)
    {
        if (z_pos[k+1] > z)
        {
            z_0 = z_pos[k];
            z_1 = z_pos[k+1];
            k_0 = k;
            break;
        }
    }

    double p_z_0[4][4];
    double p_z_1[4][4];
    for (int i=0; i<=3; i++)
    {
        for (int j=0; j<=3; j++)
        {
            p_z_0[i][j] = data[i_0+i-1][j_0+j-1][k_0];
            p_z_1[i][j] = data[i_0+i-1][j_0+j-1][k_0+1];
        }
    }

    x = (x-x_0)/(x_1-x_0);
    y = (y-y_0)/(y_1-y_0);
    z = (z-z_0)/(z_1-z_0);

    double bicubic_result_0 = bicubic_interpolate(p_z_0,x,y);
    double bicubic_result_1 = bicubic_interpolate(p_z_1,x,y);

    return bicubic_result_0+z*(bicubic_result_1-bicubic_result_0);
}


double Interpolator3D::tricubic_get_value_test (double x, double y, double z)
{
    // handling inputs outside of known range
    if (x < x_pos[0] || x > x_pos[x_pos.size()-1]) { std::cout << "Requested value outside of known range! Aborting" << std::endl; exit(0); }
    if (y < y_pos[0] || y > y_pos[y_pos.size()-1]) { std::cout << "Requested value outside of known range! Aborting" << std::endl; exit(0); }
    if (z < z_pos[0] || z > z_pos[z_pos.size()-1]) { std::cout << "Requested value outside of known range! Aborting" << std::endl; exit(0); }

    // finding corners of cuboid which x,y,z is inside of
    double x_0(0), x_1(0), y_0(0), y_1(0), z_0(0), z_1(0);
    luint i_0(0), j_0(0), k_0(0);
    for (luint i=0; i<x_pos.size()-1; i++)
    {
        if (x_pos[i+1] > x)
        {
            x_0 = x_pos[i];
            x_1 = x_pos[i+1];
            i_0 = i;
            break;
        }
    }
    for (luint j=0; j<y_pos.size()-1; j++)
    {
        if (y_pos[j+1] > y)
        {
            y_0 = y_pos[j];
            y_1 = y_pos[j+1];
            j_0 = j;
            break;
        }
    }
    for (luint k=0; k<z_pos.size()-1; k++)
    {
        if (z_pos[k+1] > z)
        {
            z_0 = z_pos[k];
            z_1 = z_pos[k+1];
            k_0 = k;
            break;
        }
    }

    double p[4][4][4];
    for (int i=0; i<=3; i++)
    {
        for (int j=0; j<=3; j++)
        {
            for (int k=0; k<=3; k++)
            {
                p[i][j][k] = data[i+i_0-1][j+j_0-1][k+k_0-1];
            }
        }
    }

    x = (x-x_0)/(x_1-x_0);
    y = (y-y_0)/(y_1-y_0);
    z = (z-z_0)/(z_1-z_0);

    double tricubic_result = tricubic_interpolate(p,x,y,z);

    return tricubic_result;
}