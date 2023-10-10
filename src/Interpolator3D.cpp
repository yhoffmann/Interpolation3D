#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cstring>
#include "../include/Interpolator3D.hpp"
#include "../external/easy-progress-monitor/include/ProgressMonitor.hpp"
#include <omp.h>


#define _INDEX(i,j,k) (i)*(n_y+3)*(n_z+3)+(j)*(n_z+3)+k


double Interpolator3D::pos_of_grid_point (Dir dir, int i, const DataGenerationConfig* config) const
{
    GridSpacing grid_spacing;
    int n;
    double min, max;
    double k;

    if (dir == Dir::x)
    {
        grid_spacing = config->x_grid_spacing;
        n = config->n_x;
        min = config->x_min;
        max = config->x_max;
        k = config->x_exp_grid_spacing_parameter;
    }
    else if (dir == Dir::y)
    {
        grid_spacing = config->y_grid_spacing;
        n = config->n_y;
        min = config->y_min;
        max = config->y_max;
        k = config->y_exp_grid_spacing_parameter;
    }
    else // if (dir == Dir::z)
    {
        grid_spacing = config->z_grid_spacing;
        n = config->n_z;
        min = config->z_min;
        max = config->z_max;
        k = config->z_exp_grid_spacing_parameter;
    }

    if (k==0.0)
        grid_spacing = Linear;

    if (grid_spacing == Linear)
        return min+(max-min)*(double(i))/(double(n-1));
    else if (grid_spacing == Exponential)
        return min+(max-min)*( exp( M_LN2*double(i)/double(n-1)*k )-1.0 )/( std::pow(2.0, k)-1.0 );
    else
        exit(-1);
}


void Interpolator3D::safe_delete_data_array()
{
    if (data_array)
        delete[] data_array;
    data_array = nullptr;
}


void Interpolator3D::prepare_data_array()
{
    safe_delete_data_array();
    data_array = new double [(n_x+3)*(n_y+3)*(n_z+3)]; // plus two because we need space for the additional element on each side
    if (!data_array)
        exit(319002);
}


void Interpolator3D::safe_delete_grid()
{
    if (x_pos)
        delete[] x_pos;
    if (y_pos)
        delete[] y_pos;
    if (z_pos)
        delete[] z_pos;

    x_pos = nullptr;
    y_pos = nullptr;
    z_pos = nullptr;
}


void Interpolator3D::set_grid (const DataGenerationConfig* config)
{
    safe_delete_grid();

    if (config)
    {
        n_x = config->n_x;
        n_y = config->n_y;
        n_z = config->n_z;
    }

    x_pos = new double [n_x+3]; // +3 because we need space for one more to the left and two more to the right
    y_pos = new double [n_y+3];
    z_pos = new double [n_z+3];

    if (!x_pos || !y_pos || !z_pos)
        exit(319001);

    if (config)
    {
        for (uint i=0; i<n_x; i++)
            x_pos[i+1] = pos_of_grid_point(Dir::x, i, config); // +1 because we do not fill the two outermost elements yet

        for (uint j=0; j<n_y; j++)
            y_pos[j+1] = pos_of_grid_point(Dir::y, j, config);

        for (uint k=0; k<n_z; k++)
            z_pos[k+1] = pos_of_grid_point(Dir::z, k, config);

        set_grid_outermost();
    }
}


void Interpolator3D::set_grid_outermost()
{
    x_pos[0] = x_pos[1]; // filling the three additional elements
    x_pos[n_x+1] = x_pos[n_x]+1.0;
    x_pos[n_x+2] = x_pos[n_x]+2.0;

    y_pos[0] = y_pos[1];
    y_pos[n_y+1] = y_pos[n_y]+1.0;
    y_pos[n_y+2] = y_pos[n_y]+2.0;

    z_pos[0] = z_pos[1];
    z_pos[n_z+1] = z_pos[n_z]+1.0;
    z_pos[n_z+2] = z_pos[n_z]+2.0;
}


void Interpolator3D::set_data_array_outermost()
{
    for (uint i=0; i<n_x+3; i++)
        for (uint j=0; j<n_y+3; j++)
            for (uint k=0; k<n_z+3; k++)
                {
                    uint i_temp = i;
                    if (i_temp==0)
                        i_temp = 1;
                    else if (i_temp>n_x+1)
                        i_temp = n_x;

                    uint j_temp = j;
                    if (j_temp==0)
                        j_temp = 1;
                    else if (j_temp>n_y+1)
                        j_temp = n_y;

                    uint k_temp = k;
                    if (k_temp==0)
                        k_temp = 1;
                    else if (k_temp>n_z)
                        k_temp = n_z;

                    data_array[_INDEX(i, j, k)] = data_array[_INDEX(i_temp, j_temp, k_temp)]; // this is doing a lot of unnecessary assignments as well but again, this would be too complicated to implement in a smarter way. also, this will typically only be run once so performance should not be an issue
                }
}


void Interpolator3D::export_data_old_format (const std::string& filepath) const
{
    // checking if file already exists
    std::ifstream file_check(filepath);
    if (file_check)
    {
        std::cout << "This file already exists and would be overwritten if data was exported. Do you want to overwrite the file? (y/n)" << std::endl;
        std::string answer;
        std::cin >> answer;
        if (answer!="y")
        {
            file_check.close();
            std::cout << "Aborting" << std::endl;
            return;
        } 
    }
    file_check.close();

    std::cout << "Exporting data_array to file..." << std::endl;

    std::ofstream out;
    out.open(filepath);
    if (!out.is_open())
    {
        std::cerr << "Could not open given file. Aborting" << std::endl;
        exit(319003);
    }

    // printing n_x, n_y, n_z into the first line
    out << "#N " << n_x << " " << n_y << " " << n_z << std::endl;

    for (uint i=0; i<n_x; i++)
    {
        for (uint j=0; j<n_y; j++)
        {
            for (uint k=0; k<n_z; k++)
            {
                out << std::setprecision(10) << i << " " << j << " " << k << " " << x_pos[i+1] << " " << y_pos[j+1] << " " << z_pos[k+1] << " " << data_array[_INDEX(i+1, j+1, k+1)] << std::endl;
            }
            out << std::endl;
        }
        out << std::endl;
    }

    out.close();

    std::cout << "Finished exporting data_array to file" << std::endl;
}


void Interpolator3D::import_data_old_format (const std::string& filepath)
{
    std::cout << "Importing data_array IN THE OLD FORMAT from file..." << std::endl;

    std::ifstream in;
    in.open(filepath);
    if (!in.is_open())
    {
        std::cerr << "Could not open given file. Aborting" << std::endl; 
        exit(319004);
    }

    std::string line;
    std::vector<std::string> line_vec;

    // getting n_x, n_y, n_z from the first line
    std::getline(in, line);
    std::istringstream ss(line);
    std::string element;
    while (std::getline(ss, element, ' '))
        line_vec.push_back(element);

    if (line_vec[0] != "#N")
    {
        std::cerr << "Wrong data_array format, use generate_data() to generate data_array in the desired format. Aborting!" << std::endl;
        exit(319005);
    }

    n_x = std::stoi(line_vec[1]);
    n_y = std::stoi(line_vec[2]);
    n_z = std::stoi(line_vec[3]);
    
    set_grid(nullptr); // allocate memory for position arrays but leave the position arrays uninitialized

    prepare_data_array();

    uint i,j,k;
    while (std::getline(in, line))
    {
        if (line != "")
        {
            std::istringstream ss(line);
            std::string element;
            line_vec.clear();
            while (std::getline(ss, element, ' '))
                line_vec.push_back(element);

            i = std::stoi(line_vec[0]);
            j = std::stoi(line_vec[1]);
            k = std::stoi(line_vec[2]);

            x_pos[i] = std::stod(line_vec[3]); // these assignments are done unnecessarily often but I cant be bothered to do this smarter
            y_pos[j] = std::stod(line_vec[4]);
            z_pos[k] = std::stod(line_vec[5]);

            data_array[_INDEX(i, j, k)] = std::stod(line_vec[6]);
        }
    }
    in.close();

    std::cout << "Finished importing data_array from file" << std::endl;
}


void Interpolator3D::export_data (const std::string& filepath) const
{
    std::ifstream file_check(filepath);
    if (file_check)
    {
        std::cout << "This file already exists and would be overwritten if data was exported. Do you want to overwrite the file? (y/n)" << std::endl;
        std::string answer;
        std::cin >> answer;
        if (answer!="y")
        {
            file_check.close();
            std::cout << "Aborting" << std::endl;
            return;
        } 
    } // TODO fix potential error that might occur when the file gets changed (for example renamed) while already opened by this function
    file_check.close();

    std::cout << "Exporting data_array to file..." << std::endl;

    std::ofstream out;
    out.open(filepath);
    if (!out.is_open())
    {
        std::cerr << "Could not open given file. Aborting" << std::endl;
        exit(319003);
    }

    out << "#n " << n_x << " " << n_y << " " << n_z << std::endl; // printing n_x, n_y, n_z into the first line

    for (uint i=0; i<n_x; i++)
        out << std::setprecision(10) << x_pos[i+1] << std::endl;
    out << std::endl; // used for detection of where x_pos ends and y_pos starts

    for (uint j=0; j<n_y; j++)
        out << std::setprecision(10) << y_pos[j+1] << std::endl;
    out << std::endl;

    for (uint k=0; k<n_z; k++)
        out << std::setprecision(10) << z_pos[k+1] << std::endl;
    out << std::endl;

    for (uint i=1; i<n_x+1; i++)
        for (uint j=1; j<n_y+1; j++)
            for (uint k=1; k<n_z+1; k++)
                out << std::setprecision(10) << data_array[_INDEX(i, j, k)] << std::endl;

    out.close();

    std::cout << "Finished exporting data_array to file" << std::endl;
}


void Interpolator3D::import_data (const std::string& filepath)
{
    std::cout << "Importing data_array from file..." << std::endl;

    std::ifstream in;
    in.open(filepath);
    if (!in.is_open())
    {
        std::cerr << "Could not open given file. Aborting" << std::endl; 
        exit(319004);
    }

    std::string line;
    std::vector<std::string> line_vec;

    // getting n_x, n_y, n_z from the first line
    std::getline(in, line);
    std::istringstream ss(line);
    std::string element;
    while (std::getline(ss, element, ' '))
        line_vec.push_back(element);

    if (line_vec[0] == "#N")
    {
        std::cerr << "=========="
        "\n"
        "You are using the new function to import data that is in the old format. This will work but you should switch to new data format by just exporting data again with export_data. Old format will remain available for export with export_data_old_format."
        "\n\nAcknowledge with enter, the program will continue as usual."
        "\n=========" << std::endl;
        std::cin.get();
        import_data_old_format(filepath);
        return;
    }
    if (line_vec[0] != "#n")
    {
        std::cerr << "Wrong data_array format, use generate_data() to generate data_array in the desired format. Aborting!" << std::endl;
        exit(319005);
    }

    n_x = std::stoi(line_vec[1]);
    n_y = std::stoi(line_vec[2]);
    n_z = std::stoi(line_vec[3]);
    
    set_grid(nullptr); // allocate memory for position arrays but leave the position arrays uninitialized

    prepare_data_array();

    uint i(0), j(0), k(0);
    while (std::getline(in, line) && line != "")
        x_pos[++i] = std::stod(line);

    while (std::getline(in, line) && line != "")
        y_pos[++j] = std::stod(line);
    
    while (std::getline(in, line) && line != "")
        z_pos[++k] = std::stod(line);

    set_grid_outermost();

    for (i=1; i<n_x+1; i++)
        for (j=1; j<n_y+1; j++)
            for (k=1; k<n_z+1; k++)
            {
                std::getline(in, line);
                data_array[_INDEX(i, j, k)] = std::stod(line);
            }

    set_data_array_outermost();

    in.close();

    std::cout << "Finished importing data_array from file" << std::endl;
}


void Interpolator3D::generate_data (std::function<double (double,double,double)> func, const DataGenerationConfig* config, bool progress_monitor)
{
    set_grid(config);
    prepare_data_array();

    ProgressMonitor pm(n_x);

    // filling data_array with function values
    #pragma omp parallel for ordered
    for (uint i=1; i<n_x+1; i++)
    {
        for (uint j=1; j<n_y+1; j++)
        {
            for (uint k=1; k<n_z+1; k++)
            {
                data_array[_INDEX(i, j, k)] = func(x_pos[i], y_pos[j], z_pos[k]);
            }
        }
        if (progress_monitor)
        {
            pm.add_finished();
            pm.print_progress_percentage();
        }
    }

    set_data_array_outermost();
}


void Interpolator3D::find_closest_lower_data_point(int& i_0, int& j_0, int& k_0, double& x, double& y, double& z) const
{
    i_0 = find_index(x_pos+1, n_x, x);
    j_0 = find_index(y_pos+1, n_y, y);
    k_0 = find_index(z_pos+1, n_z, z);

    ++i_0; // to account for the outermost values
    ++j_0;
    ++k_0;

    if (x<x_pos[1])
        x = x_pos[1];
    else if (x>x_pos[n_x])
        x = x_pos[n_x];

    if (y<y_pos[1])
        y = y_pos[1];
    else if (y>y_pos[n_y])
        y = y_pos[n_y];

    if (z<z_pos[1])
        z = z_pos[1];
    else if (z>z_pos[n_z])
        z = z_pos[n_z];
}


uint Interpolator3D::find_index (double* arr, int size, double x)
{
    int low = 0;
    int high = size-1;

    while (low <= high)
    {
        int mid = low + (high-low)/2;

        if (arr[mid] <= x)
        {
            if (mid==size-1 || arr[mid+1]>x)
                return mid;
            else
                low = mid+1;
        }
        else
            high = mid-1;
    }

    return 0;
}


double Interpolator3D::safe_get_x_pos (int i) const
{
    if (i<0)
        return x_pos[0]+double(i);
    else if (i>int(n_x-1))
        return x_pos[n_x-1]+double(i-int(n_x-1));
    else
        return x_pos[i];
}


double Interpolator3D::safe_get_y_pos (int i) const
{
    if (i<0)
        return y_pos[0]+double(i);
    else if (i>int(n_y-1))
        return y_pos[n_y-1]+double(i-int(n_y-1));
    else
        return y_pos[i];
}


double Interpolator3D::safe_get_z_pos (int i) const
{
    if (i<0)
        return z_pos[0]+double(i);
    else if (i>int(n_z-1))
        return z_pos[n_z-1]+double(i-int(n_z-1));
    else
        return z_pos[i];
}


double Interpolator3D::safe_get_data_point (int i, int j, int k) const
{
    if (i<0)
        i = 0;
    else if (i>int(n_x-1)) 
        i = n_x-1;

    if (j<0)
        j = 0;
    else if (j>int(n_y-1))
        j = n_y-1;

    if (k<0)
        k = 0;
    else if (k>int(n_z-1))
        k = n_z-1;

    return data_array[_INDEX(i, j, k)];
}


double Interpolator3D::unicubic_interpolate (double p[4], double t[4], double z)
{
    double p2_m_p0_div_t2_m_t0 = (p[2]-p[0])/(1.0-t[0]);
    double p3_m_p1_div_t3_m_t1 = (p[3]-p[1])/t[3];
    double p1_m_p2 = p[1]-p[2];

    return p[1] + z*p2_m_p0_div_t2_m_t0 + z*z*(-3.0*p1_m_p2-2.0*p2_m_p0_div_t2_m_t0-p3_m_p1_div_t3_m_t1) + z*z*z*(2.0*p1_m_p2+p2_m_p0_div_t2_m_t0+p3_m_p1_div_t3_m_t1);
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


double Interpolator3D::get_interp_value_bicubic_unilinear (double x, double y, double z) const
{
    int i_0, j_0, k_0;

    find_closest_lower_data_point(i_0, j_0, k_0, x, y, z);

    double p_z_0[4][4];
    double p_z_1[4][4];
    for (int i=0; i<4; i++)
    {
        for (int j=0; j<4; j++)
        {
            p_z_0[i][j] = data_array[_INDEX(i+i_0-1, j+j_0-1, k_0)];
            p_z_1[i][j] = data_array[_INDEX(i+i_0-1, j+j_0-1, k_0+1)];
        }
    }

    double t_x[4];
    double t_y[4];
    for (int i=0; i<4; i++)
    {
        t_x[i] = x_pos[i+i_0-1];
        t_y[i] = y_pos[i+j_0-1];
    }

    double one_div_t_x_2_t_x_1 = 1.0/(t_x[2]-t_x[1]); // this looks like this because of performance advantages
    double one_div_t_y_2_t_y_1 = 1.0/(t_y[2]-t_y[1]);

    x = (x-t_x[1])*one_div_t_x_2_t_x_1;
    y = (y-t_y[1])*one_div_t_y_2_t_y_1;
    z = (z-z_pos[k_0])/(z_pos[k_0+1]-z_pos[k_0]);

    t_x[0] = (t_x[0]-t_x[1])*one_div_t_x_2_t_x_1;
    t_x[3] = (t_x[3]-t_x[1])*one_div_t_x_2_t_x_1;
    t_x[1] = 0.0;
    t_x[2] = 1.0;

    t_y[0] = (t_y[0]-t_y[1])*one_div_t_y_2_t_y_1;
    t_y[3] = (t_y[3]-t_y[1])*one_div_t_y_2_t_y_1;
    t_y[1] = 0.0;
    t_y[2] = 1.0;

    double bicbuic_result_0 = bicubic_interpolate(p_z_0, t_x, t_y, x, y);
    double bicbuic_result_1 = bicubic_interpolate(p_z_1, t_x, t_y, x, y);

    double bicubic_unilinear_result = bicbuic_result_0+z*(bicbuic_result_1-bicbuic_result_0);

    return bicubic_unilinear_result;
}


double Interpolator3D::get_interp_value_tricubic (double x, double y, double z) const
{
    int i_0, j_0, k_0;

    find_closest_lower_data_point(i_0, j_0, k_0, x, y, z);

    double p[4][4][4];

    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            for (int k=0; k<4; k++)
                p[i][j][k] = data_array[_INDEX(i+i_0-1, j+j_0-1, k+k_0-1)];

    double t_x[4];
    double t_y[4];
    double t_z[4];

    for (int i=0; i<4; i++)
    {
        t_x[i] = x_pos[i+i_0-1];
        t_y[i] = y_pos[i+j_0-1];
        t_z[i] = z_pos[i+k_0-1];
    }

    double one_div_t_x_2_t_x_1 = 1.0/(t_x[2]-t_x[1]); // this looks like this because of performance advantages
    double one_div_t_y_2_t_y_1 = 1.0/(t_y[2]-t_y[1]);
    double one_div_t_z_2_t_z_1 = 1.0/(t_z[2]-t_z[1]);

    x = (x-t_x[1])*one_div_t_x_2_t_x_1;
    y = (y-t_y[1])*one_div_t_y_2_t_y_1;
    z = (z-t_z[1])*one_div_t_z_2_t_z_1;

    t_x[0] = (t_x[0]-t_x[1])*one_div_t_x_2_t_x_1;
    t_x[3] = (t_x[3]-t_x[1])*one_div_t_x_2_t_x_1;
    t_x[1] = 0.0;
    t_x[2] = 1.0;

    t_y[0] = (t_y[0]-t_y[1])*one_div_t_y_2_t_y_1;
    t_y[3] = (t_y[3]-t_y[1])*one_div_t_y_2_t_y_1;
    t_y[1] = 0.0;
    t_y[2] = 1.0;

    t_z[0] = (t_z[0]-t_z[1])*one_div_t_z_2_t_z_1;
    t_z[3] = (t_z[3]-t_z[1])*one_div_t_z_2_t_z_1;
    t_z[1] = 0.0;
    t_z[2] = 1.0;

    double tricubic_result = tricubic_interpolate(p, t_x, t_y, t_z, x, y, z);

    return tricubic_result;
}


Interpolator3D::Interpolator3D()
{
}


Interpolator3D::Interpolator3D (const std::string& filepath)
{
    import_data(filepath);
}


Interpolator3D::Interpolator3D (const Interpolator3D& other)
{
    n_x = other.n_x;
    n_y = other.n_y;
    n_z = other.n_z;

    prepare_data_array();
    set_grid(nullptr);

    std::copy(other.data_array, other.data_array+(n_x+3)*(n_y+3)*(n_z+3), data_array);
    std::copy(other.x_pos, other.x_pos+n_x+3, x_pos);
    std::copy(other.y_pos, other.y_pos+n_y+3, y_pos);
    std::copy(other.z_pos, other.z_pos+n_z+3, z_pos);
}


Interpolator3D& Interpolator3D::operator= (const Interpolator3D& other)
{
    if (this==&other)
        return *this;

    n_x = other.n_x;
    n_y = other.n_y;
    n_z = other.n_z;

    prepare_data_array();
    set_grid(nullptr);

    std::copy(other.data_array, other.data_array+(n_x+3)*(n_y+3)*(n_z+3), data_array);
    std::copy(other.x_pos, other.x_pos+n_x+3, x_pos);
    std::copy(other.y_pos, other.y_pos+n_y+3, y_pos);
    std::copy(other.z_pos, other.z_pos+n_z+3, z_pos);

    return *this;
}


Interpolator3D& Interpolator3D::operator= (Interpolator3D&& other)
{
    if (this==&other)
        return *this;

    memcpy(this, &other, sizeof(Interpolator3D));

    other.data_array = nullptr;
    other.x_pos = nullptr;
    other.y_pos = nullptr;
    other.z_pos = nullptr;

    return *this;
}


Interpolator3D::~Interpolator3D()
{
    safe_delete_grid();
    
    safe_delete_data_array();
}