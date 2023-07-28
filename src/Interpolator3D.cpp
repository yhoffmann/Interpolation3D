#include "../include/Interpolator3D.h"
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>


double Interpolator3D::inputs_for_pos(XYZ xyz, luint index, DataGenerationConfig& config)
{
    std::string grid_spacing;
    luint n;
    double min, max;

    if (xyz == x) { grid_spacing = config.x_grid_spacing; n = config.n_x; min = config.x_min; max = config.x_max; }
    else if (xyz == y) { grid_spacing = config.y_grid_spacing; n = config.n_y; min = config.y_min; max = config.y_max; }
    else if (xyz == z) { grid_spacing = config.z_grid_spacing; n = config.n_z; min = config.z_min; max = config.z_max; }

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
        x_pos[i] = inputs_for_pos(x,i,config);
    }
    for (luint j=0; j<n_y; j++)
    {
        y_pos[j] = inputs_for_pos(y,j,config);
    }
    for (luint k=0; k<n_z; k++)
    {
        z_pos[k] = inputs_for_pos(z,k,config);
    }
}


void Interpolator3D::print_data_to_file (std::string filepath) 
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
            for (int k=0; k<data[0][0].size(); k++)
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
        if (line != "\n")
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

            x_pos[i] = std::stod(line_vec[3]); // these assignments are done unnessecarily often but I cant be bothered to do this smarter
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
    double x_0, x_1, y_0, y_1, z_0, z_1;
    luint i_0, j_0, k_0;
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