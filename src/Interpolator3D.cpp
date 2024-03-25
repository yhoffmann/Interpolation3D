#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <omp.h>
#include <algorithm>
#include "../include/Interpolator3D.hpp"
#include "../external/easy-progress-monitor/include/ProgressMonitor.hpp"


#define _INDEX(i, j, k) ((i)*(m_ny+3)*(m_nz+3)+(j)*(m_nz+3)+k)


double Interpolator3D::pos_of_grid_point (Dir dir, int i, const DataGenerationConfig* config) const
{
    GridSpacing grid_spacing;
    int n;
    double min, max;
    double k;

    if (dir == Dir::x)
    {
        grid_spacing = config->x_grid_spacing;
        n = config->nx;
        min = config->x_min;
        max = config->x_max;
        k = config->x_exp_grid_spacing_parameter;
    }
    else if (dir == Dir::y)
    {
        grid_spacing = config->y_grid_spacing;
        n = config->ny;
        min = config->y_min;
        max = config->y_max;
        k = config->y_exp_grid_spacing_parameter;
    }
    else // if (dir == Dir::z)
    {
        grid_spacing = config->z_grid_spacing;
        n = config->nz;
        min = config->z_min;
        max = config->z_max;
        k = config->z_exp_grid_spacing_parameter;
    }

    if (k==0.0)
        grid_spacing = Linear;

    if (grid_spacing == Linear)
        return min+(max-min)*(double(i))/(double(n-1));
    else // if (grid_spacing == Exponential)
        return min+(max-min)*( exp( M_LN2*double(i)/double(n-1)*k )-1.0 )/( std::pow(2.0, k)-1.0 );
}


void Interpolator3D::safe_delete_m_data()
{
    if (m_data)
    {
        delete[] m_data;
        m_data = nullptr;
    }
}


void Interpolator3D::prepare_m_data()
{
    safe_delete_m_data();
    m_data = new(std::nothrow) double [(m_nx+3)*(m_ny+3)*(m_nz+3)]; // plus two because we need space for the additional element on each side
    if (!m_data)
        exit(50);
}


void Interpolator3D::safe_delete_grid()
{
    if (m_x)
    {
        delete[] m_x;
        m_x = nullptr;
    }
    if (m_y)
    {
        delete[] m_y;
        m_y = nullptr;
    }
    if (m_z)
    {
        delete[] m_z;
        m_z = nullptr;
    }
}


void Interpolator3D::set_grid (const DataGenerationConfig* config)
{
    safe_delete_grid();

    if (config)
    {
        m_nx = config->nx;
        m_ny = config->ny;
        m_nz = config->nz;
    }

    m_x = new(std::nothrow) double [m_nx+3]; // +3 because we need space for one more to the left and two more to the right
    m_y = new(std::nothrow) double [m_ny+3];
    m_z = new(std::nothrow) double [m_nz+3];

    if (!m_x || !m_y || !m_z)
        exit(51);

    if (config)
    {
        for (uint i=0; i<m_nx; i++)
            m_x[i+1] = pos_of_grid_point(Dir::x, i, config); // +1 because we do not fill the two outermost elements yet

        for (uint j=0; j<m_ny; j++)
            m_y[j+1] = pos_of_grid_point(Dir::y, j, config);

        for (uint k=0; k<m_nz; k++)
            m_z[k+1] = pos_of_grid_point(Dir::z, k, config);

        set_grid_outermost();
    }
}


void Interpolator3D::set_grid_outermost()
{
    m_x[0] = m_x[1]-1.0; // filling the three additional elements
    m_x[m_nx+1] = m_x[m_nx]+1.0;
    m_x[m_nx+2] = m_x[m_nx+1]+1.0;

    m_y[0] = m_y[1]-1.0;
    m_y[m_ny+1] = m_y[m_ny]+1.0;
    m_y[m_ny+2] = m_y[m_ny+1]+1.0;

    m_z[0] = m_z[1]-1.0;
    m_z[m_nz+1] = m_z[m_nz]+1.0;
    m_z[m_nz+2] = m_z[m_nz+1]+1.0;
}


void Interpolator3D::set_m_data_outermost()
{
    for (uint i=0; i<m_nx+3; i++)
        for (uint j=0; j<m_ny+3; j++)
            for (uint k=0; k<m_nz+3; k++)
                {
                    uint i_temp = i;
                    if (i_temp==0)
                        i_temp = 1;
                    else if (i_temp>m_nx)
                        i_temp = m_nx;

                    uint j_temp = j;
                    if (j_temp==0)
                        j_temp = 1;
                    else if (j_temp>m_ny)
                        j_temp = m_ny;

                    uint k_temp = k;
                    if (k_temp==0)
                        k_temp = 1;
                    else if (k_temp>m_nz)
                        k_temp = m_nz;

                    m_data[_INDEX(i, j, k)] = m_data[_INDEX(i_temp, j_temp, k_temp)]; // this is doing a lot of unnecessary assignments as well but again, this would be too complicated to implement in a smarter way. also, this will typically only be run once so performance should not be an issue
                }
}


void Interpolator3D::export_data_old_format (const std::string& filepath) const
{
#ifndef _QUIET
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

    std::cout << "Exporting m_data to file..." << std::endl;
#endif

    std::ofstream out;
    out.open(filepath);
    if (!out.is_open())
    {
#ifndef _QUIET
        std::cerr << "Could not open given file. Aborting" << std::endl;
#endif
        exit(52);
    }

    // printing m_nx, m_ny, m_nz into the first line
    out << "#N " << m_nx << " " << m_ny << " " << m_nz << std::endl;

    for (uint i=0; i<m_nx; i++)
    {
        for (uint j=0; j<m_ny; j++)
        {
            for (uint k=0; k<m_nz; k++)
            {
                out << std::setprecision(10) << i << " " << j << " " << k << " " << m_x[i+1] << " " << m_y[j+1] << " " << m_z[k+1] << " " << m_data[_INDEX(i+1, j+1, k+1)] << std::endl;
            }
            out << std::endl;
        }
        out << std::endl;
    }

    out.close();
#ifndef _QUIET
    std::cout << "Finished exporting m_data to file" << std::endl;
#endif
}


void Interpolator3D::import_data_old_format (const std::string& filepath)
{
#ifndef _QUIET
    std::cout << "Importing m_data IN THE OLD FORMAT from file..." << std::endl;
#endif
    std::ifstream in;
    in.open(filepath);
    if (!in.is_open())
    {
#ifndef _QUIET
        std::cerr << "Could not open given file. Aborting" << std::endl; 
#endif
        exit(53);
    }

    std::string line;
    std::vector<std::string> line_vec;

    // getting m_nx, m_ny, m_nz from the first line
    std::getline(in, line);
    std::istringstream ss(line);
    std::string element;
    while (std::getline(ss, element, ' '))
        line_vec.push_back(element);

    if (line_vec[0] != "#N")
    {
#ifndef _QUIET
        std::cerr << "Wrong m_data format, use generate_data() to generate m_data in the desired format. Aborting!" << std::endl;
#endif
        exit(54);
    }

    m_nx = std::stoi(line_vec[1]);
    m_ny = std::stoi(line_vec[2]);
    m_nz = std::stoi(line_vec[3]);
    
    set_grid(nullptr); // allocate memory for position arrays but leave the position arrays uninitialized

    prepare_m_data();

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

            m_x[i+1] = std::stod(line_vec[3]); // these assignments are done unnecessarily often but I cant be bothered to do this smarter
            m_y[j+1] = std::stod(line_vec[4]);
            m_z[k+1] = std::stod(line_vec[5]);

            m_data[_INDEX(i+1, j+1, k+1)] = std::stod(line_vec[6]);
        }
    }

    in.close();

    set_grid_outermost();
    set_m_data_outermost();
#ifndef _QUIET
    std::cout << "Finished importing m_data from file" << std::endl;
#endif
}


void Interpolator3D::export_data (const std::string& filepath) const
{
#ifndef _QUIET
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

    std::cout << "Exporting m_data to file..." << std::endl;
#endif
    std::ofstream out;
    out.open(filepath);
    if (!out.is_open())
    {
#ifndef _QUIET
        std::cerr << "Could not open given file. Aborting" << std::endl;
#endif
        exit(52);
    }

    out << "#n " << m_nx << " " << m_ny << " " << m_nz << std::endl; // printing m_nx, m_ny, m_nz into the first line

    for (uint i=0; i<m_nx; i++)
        out << std::setprecision(10) << m_x[i+1] << std::endl;
    out << std::endl; // used for detection of where m_x ends and m_y starts

    for (uint j=0; j<m_ny; j++)
        out << std::setprecision(10) << m_y[j+1] << std::endl;
    out << std::endl;

    for (uint k=0; k<m_nz; k++)
        out << std::setprecision(10) << m_z[k+1] << std::endl;
    out << std::endl;

    for (uint i=1; i<m_nx+1; i++)
        for (uint j=1; j<m_ny+1; j++)
            for (uint k=1; k<m_nz+1; k++)
                out << std::setprecision(10) << m_data[_INDEX(i, j, k)] << std::endl;

    out.close();
#ifndef _QUIET
    std::cout << "Finished exporting m_data to file" << std::endl;
#endif
}


void Interpolator3D::import_data (const std::string& filepath)
{
#ifndef _QUIET
    std::cout << "Importing m_data from file..." << std::endl;
#endif

    std::ifstream in;
    in.open(filepath);
    if (!in.is_open())
    {
#ifndef _QUIET
        std::cerr << "Could not open given file. Aborting" << std::endl;
#endif
        exit(53);
    }

    std::string line;
    std::vector<std::string> line_vec;

    // getting m_nx, m_ny, m_nz from the first line
    std::getline(in, line);
    std::istringstream ss(line);
    std::string element;
    while (std::getline(ss, element, ' '))
        line_vec.push_back(element);

    if (line_vec[0] == "#N")
    {
#ifndef _QUIET
        std::cerr << "=========="
        "\n"
        "You are using the new function to import data that is in the old format. This will work but you should switch to new data format by just exporting data again with export_data. Old format will remain available for export with export_data_old_format."
        "\n\nAcknowledge with enter, the program will continue as usual."
        "\n=========" << std::endl;
        std::cin.get();
#endif
        import_data_old_format(filepath);
        return;
    }
    if (line_vec[0] != "#n")
    {
#ifndef _QUIET
        std::cerr << "Wrong m_data format, use generate_data() to generate m_data in the desired format. Aborting!" << std::endl;
#endif
        exit(54);
    }

    m_nx = std::stoi(line_vec[1]);
    m_ny = std::stoi(line_vec[2]);
    m_nz = std::stoi(line_vec[3]);
    
    set_grid(nullptr); // allocate memory for position arrays but leave the position arrays uninitialized

    prepare_m_data();

    uint i(0), j(0), k(0);
    while (std::getline(in, line) && line != "")
        m_x[++i] = std::stod(line);

    while (std::getline(in, line) && line != "")
        m_y[++j] = std::stod(line);
    
    while (std::getline(in, line) && line != "")
        m_z[++k] = std::stod(line);

    set_grid_outermost();

    for (i=1; i<m_nx+1; i++)
        for (j=1; j<m_ny+1; j++)
            for (k=1; k<m_nz+1; k++)
            {
                std::getline(in, line);
                m_data[_INDEX(i, j, k)] = std::stod(line);
            }

    in.close();

    set_m_data_outermost();
#ifndef _QUIET
    std::cout << "Finished importing m_data from file" << std::endl;
#endif
}


void Interpolator3D::generate_data (std::function<double (double,double,double)> func, const DataGenerationConfig* config, bool progress_monitor)
{
    set_grid(config);
    prepare_m_data();

    ProgressMonitor pm(m_nx);

    // filling m_data with function values
    #pragma omp parallel for ordered
    for (uint i=1; i<m_nx+1; i++)
    {
        for (uint j=1; j<m_ny+1; j++)
        {
            for (uint k=1; k<m_nz+1; k++)
            {
                m_data[_INDEX(i, j, k)] = func(m_x[i], m_y[j], m_z[k]);
            }
        }
#ifndef _QUIET
        if (progress_monitor)
        {
            pm.add_finished();
            pm.print_progress_percentage();
        }
#endif
    }

    set_m_data_outermost();
}


void Interpolator3D::find_closest_lower_data_point(int& i_0, int& j_0, int& k_0, double& x, double& y, double& z) const
{
    i_0 = find_index(m_x+1, m_nx, x);
    j_0 = find_index(m_y+1, m_ny, y);
    k_0 = find_index(m_z+1, m_nz, z);

    // ++i_0; // +1 removed here and added everywhere in interpolation functions
    // ++j_0;
    // ++k_0;

    if (x<m_x[1])
        x = m_x[1];
    else if (x>m_x[m_nx])
        x = m_x[m_nx];

    if (y<m_y[1])
        y = m_y[1];
    else if (y>m_y[m_ny])
        y = m_y[m_ny];

    if (z<m_z[1])
        z = m_z[1];
    else if (z>m_z[m_nz])
        z = m_z[m_nz];
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
            p_z_0[i][j] = m_data[_INDEX(i+i_0, j+j_0, k_0+1)];
            p_z_1[i][j] = m_data[_INDEX(i+i_0, j+j_0, k_0+2)];
        }
    }

    double t_x[4];
    double t_y[4];
    for (int i=0; i<4; i++)
        t_x[i] = m_x[i+i_0];
    for (int i=0; i<4; i++)
        t_y[i] = m_y[i+j_0];

    double one_div_t_x_2_t_x_1 = 1.0/(t_x[2]-t_x[1]); // this looks like this because of performance advantages
    double one_div_t_y_2_t_y_1 = 1.0/(t_y[2]-t_y[1]);

    x = (x-t_x[1])*one_div_t_x_2_t_x_1;
    y = (y-t_y[1])*one_div_t_y_2_t_y_1;
    z = (z-m_z[k_0])/(m_z[k_0+1]-m_z[k_0]);

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

    return bicbuic_result_0+z*(bicbuic_result_1-bicbuic_result_0);
}


double Interpolator3D::get_interp_value_tricubic (double x, double y, double z) const
{
    int i_0, j_0, k_0;

    find_closest_lower_data_point(i_0, j_0, k_0, x, y, z);

    double p[4][4][4];

    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            for (int k=0; k<4; k++)
                p[i][j][k] = m_data[_INDEX(i+i_0, j+j_0, k+k_0)];

    double t_x[4];
    double t_y[4];
    double t_z[4];

    for (int i=0; i<4; i++)
        t_x[i] = m_x[i+i_0];
    for (int i=0; i<4; i++)
        t_y[i] = m_y[i+j_0];
    for (int i=0; i<4; i++)
        t_z[i] = m_z[i+k_0];

    double one_div_t_x_2_m_t_x_1 = 1.0/(t_x[2]-t_x[1]); // this looks like this because of performance advantages
    double one_div_t_y_2_m_t_y_1 = 1.0/(t_y[2]-t_y[1]);
    double one_div_t_z_2_m_t_z_1 = 1.0/(t_z[2]-t_z[1]);

    x = (x-t_x[1])*one_div_t_x_2_m_t_x_1;
    y = (y-t_y[1])*one_div_t_y_2_m_t_y_1;
    z = (z-t_z[1])*one_div_t_z_2_m_t_z_1;

    t_x[0] = (t_x[0]-t_x[1])*one_div_t_x_2_m_t_x_1;
    t_x[3] = (t_x[3]-t_x[1])*one_div_t_x_2_m_t_x_1;
    t_x[1] = 0.0;
    t_x[2] = 1.0;

    t_y[0] = (t_y[0]-t_y[1])*one_div_t_y_2_m_t_y_1;
    t_y[3] = (t_y[3]-t_y[1])*one_div_t_y_2_m_t_y_1;
    t_y[1] = 0.0;
    t_y[2] = 1.0;

    t_z[0] = (t_z[0]-t_z[1])*one_div_t_z_2_m_t_z_1;
    t_z[3] = (t_z[3]-t_z[1])*one_div_t_z_2_m_t_z_1;
    t_z[1] = 0.0;
    t_z[2] = 1.0;

    return tricubic_interpolate(p, t_x, t_y, t_z, x, y, z);
}


Interpolator3D::Interpolator3D()
{
}


Interpolator3D::Interpolator3D (const std::string& filepath)
{
    import_data(filepath);
}


Interpolator3D::Interpolator3D (const Interpolator3D& other)
    : m_nx(other.m_nx)
    , m_ny(other.m_ny)
    , m_nz(other.m_nz)
{

    prepare_m_data();
    set_grid(nullptr);

    std::copy(other.m_data, other.m_data+(m_nx+3)*(m_ny+3)*(m_nz+3), m_data);
    std::copy(other.m_x, other.m_x+m_nx+3, m_x);
    std::copy(other.m_y, other.m_y+m_ny+3, m_y);
    std::copy(other.m_z, other.m_z+m_nz+3, m_z);
}


Interpolator3D::Interpolator3D (Interpolator3D&& other)
    : m_nx(other.m_nx)
    , m_ny(other.m_ny)
    , m_nz(other.m_nz)
    , m_data(other.m_data)
    , m_x(other.m_x)
    , m_y(other.m_y)
    , m_z(other.m_z)
{
    other.m_data = nullptr;
    other.m_x = nullptr;
    other.m_y = nullptr;
    other.m_z = nullptr;
}


Interpolator3D& Interpolator3D::operator= (const Interpolator3D& other)
{
    if (this==&other)
        return *this;

    m_nx = other.m_nx;
    m_ny = other.m_ny;
    m_nz = other.m_nz;

    prepare_m_data();
    set_grid(nullptr);

    std::copy(other.m_data, other.m_data+(m_nx+3)*(m_ny+3)*(m_nz+3), m_data);
    std::copy(other.m_x, other.m_x+m_nx+3, m_x);
    std::copy(other.m_y, other.m_y+m_ny+3, m_y);
    std::copy(other.m_z, other.m_z+m_nz+3, m_z);

    return *this;
}


Interpolator3D& Interpolator3D::operator= (Interpolator3D&& other)
{
    if (this==&other)
        return *this;

    std::copy(&other, &other+1, this);

    other.m_data = nullptr;
    other.m_x = nullptr;
    other.m_y = nullptr;
    other.m_z = nullptr;

    return *this;
}


Interpolator3D::~Interpolator3D()
{
    safe_delete_grid();
    
    safe_delete_m_data();
}


#undef _INDEX