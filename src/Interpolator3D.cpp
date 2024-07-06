#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <gsl/gsl_math.h>
#include "../include/Interpolator3D.hpp"
#include "../external/easy-progress-monitor/include/ProgressMonitor.hpp"
#include "../include/CoeffMatrix.hpp"

#define _INDEX(i, j, k) ((i)*(m_ny+3)*(m_nz+3)+(j)*(m_nz+3)+k)

typedef int InterpFormatCheck;
const InterpFormatCheck FORMAT_CHECK_NUM = 1906853104;


double Interpolator3D::pos_of_grid_point (Dir dir, int i, const DataGenerationConfig* config) const
{
    GridSpacing grid_spacing;
    int n;
    double min, max;
    double k;

    switch (dir)
    {
        case Dir::X:
            grid_spacing = config->x_grid_spacing;
            n = config->nx;
            min = config->x_min;
            max = config->x_max;
            k = config->x_exp_grid_spacing_parameter;
        break;

        case Dir::Y:
            grid_spacing = config->y_grid_spacing;
            n = config->ny;
            min = config->y_min;
            max = config->y_max;
            k = config->y_exp_grid_spacing_parameter;
        break;

        case Dir::Z:
        default:
            grid_spacing = config->z_grid_spacing;
            n = config->nz;
            min = config->z_min;
            max = config->z_max;
            k = config->z_exp_grid_spacing_parameter;
        break;
    }

    if (k==0.0)
        grid_spacing = Linear;

    if (grid_spacing == Linear)
        return min+(max-min)*double(i)/double(n-1);
    else // if (grid_spacing == Exponential)
        return min+(max-min)*( exp( M_LN2*double(i)/double(n-1)*k )-1.0 )/( std::pow(2.0, k)-1.0 );
}


void Interpolator3D::safe_delete_data()
{
    if (m_data!=nullptr)
        delete[] m_data;
}


void Interpolator3D::prepare_data()
{
    safe_delete_data();
    m_data = new(std::nothrow) double [(m_nx+3)*(m_ny+3)*(m_nz+3)];
    if (!m_data)
        exit(50);
}


void Interpolator3D::safe_delete_cached_coeffs()
{
    if (m_cached_coeffs!=nullptr)
        delete[] m_cached_coeffs;
}


void Interpolator3D::prepare_cached_coeffs()
{
    safe_delete_cached_coeffs();
    m_cached_coeffs = new(std::nothrow) Coeffs[m_nx*m_ny*m_nz];

    if (!m_cached_coeffs)
        exit(56);
}


void Interpolator3D::cache_coeffs()
{
    prepare_cached_coeffs();

    for (uint i=0; i<m_nx; i++)
        for (uint j=0; j<m_ny; j++)
            for (uint k=0; k<m_nz; k++)
                set_single_cell_coeffs(i, j, k);
}


void Interpolator3D::set_single_cell_coeffs(uint i_0, uint j_0, uint k_0)
{
    double p[4][4][4];

    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            for (int k=0; k<4; k++)
                p[i][j][k] = m_data[_INDEX(i+i_0, j+j_0, k+k_0)];

    double tx[4];
    double ty[4];
    double tz[4];
    for (uint i=0; i<4; i++)
    {
        tx[i] = m_x[i_0+i];
        ty[i] = m_y[j_0+i];
        tz[i] = m_z[k_0+i];
    }

    tx[0] = (tx[0]-tx[1])/(tx[2]-tx[1]);
    tx[3] = (tx[3]-tx[1])/(tx[2]-tx[1]);
    double one_div_dx2 = 1.0/(1.0-tx[0]);
    double one_div_dx3 = 1.0/tx[3];

    ty[0] = (ty[0]-ty[1])/(ty[2]-ty[1]);
    ty[3] = (ty[3]-ty[1])/(ty[2]-ty[1]);
    double one_div_dy2 = 1.0/(1.0-ty[0]);
    double one_div_dy3 = 1.0/ty[3];

    tz[0] = (tz[0]-tz[1])/(tz[2]-tz[1]);
    tz[3] = (tz[3]-tz[1])/(tz[2]-tz[1]);
    double one_div_dz2 = 1.0/(1.0-tz[0]);
    double one_div_dz3 = 1.0/tz[3];

    // coeffs = COEFF_MATRIX * b
    Coeffs b;

    // values
    b[0] = p[1][1][1];
    b[1] = p[2][1][1];
    b[2] = p[1][2][1];
    b[3] = p[2][2][1];
    b[4] = p[1][1][2];
    b[5] = p[2][1][2];
    b[6] = p[1][2][2];
    b[7] = p[2][2][2];

    // partial x derivatives
    b[ 8] = (p[2][1][1]-p[0][1][1])*one_div_dx2;
    b[ 9] = (p[3][1][1]-p[1][1][1])*one_div_dx3;
    b[10] = (p[2][2][1]-p[0][2][1])*one_div_dx2;
    b[11] = (p[3][2][1]-p[1][2][1])*one_div_dx3;
    b[12] = (p[2][1][2]-p[0][1][2])*one_div_dx2;
    b[13] = (p[3][1][2]-p[1][1][2])*one_div_dx3;
    b[14] = (p[2][2][2]-p[0][2][2])*one_div_dx2;
    b[15] = (p[3][2][2]-p[1][2][2])*one_div_dx3;

    // partial y derivatives
    b[16] = (p[1][2][1]-p[1][0][1])*one_div_dy2;
    b[17] = (p[2][2][1]-p[2][0][1])*one_div_dy2;
    b[18] = (p[1][3][1]-p[1][1][1])*one_div_dy3;
    b[19] = (p[2][3][1]-p[2][1][1])*one_div_dy3;
    b[20] = (p[1][2][2]-p[1][0][2])*one_div_dy2;
    b[21] = (p[2][2][2]-p[2][0][2])*one_div_dy2;
    b[22] = (p[1][3][2]-p[1][1][2])*one_div_dy3;
    b[23] = (p[2][3][2]-p[2][1][2])*one_div_dy3;

    // partial z derivatives
    b[24] = (p[1][1][2]-p[1][1][0])*one_div_dz2;
    b[25] = (p[2][1][2]-p[2][1][0])*one_div_dz2;
    b[26] = (p[1][2][2]-p[1][2][0])*one_div_dz2;
    b[27] = (p[2][2][2]-p[2][2][0])*one_div_dz2;
    b[28] = (p[1][1][3]-p[1][1][1])*one_div_dz3;
    b[29] = (p[2][1][3]-p[2][1][1])*one_div_dz3;
    b[30] = (p[1][2][3]-p[1][2][1])*one_div_dz3;
    b[31] = (p[2][2][3]-p[2][2][1])*one_div_dz3;

    // partial xy derivatives
    b[32] = (p[2][2][1]-p[2][0][1]-(p[0][2][1]-p[0][0][1]))*one_div_dx2*one_div_dy2;
    b[33] = (p[3][2][1]-p[3][0][1]-(p[1][2][1]-p[1][0][1]))*one_div_dx3*one_div_dy2;
    b[34] = (p[2][3][1]-p[2][1][1]-(p[0][3][1]-p[0][1][1]))*one_div_dx2*one_div_dy3;
    b[35] = (p[3][3][1]-p[3][1][1]-(p[1][3][1]-p[1][1][1]))*one_div_dx3*one_div_dy3;
    b[36] = (p[2][2][2]-p[2][0][2]-(p[0][2][2]-p[0][0][2]))*one_div_dx2*one_div_dy2;
    b[37] = (p[3][2][2]-p[3][0][2]-(p[1][2][2]-p[1][0][2]))*one_div_dx3*one_div_dy2;
    b[38] = (p[2][3][2]-p[2][1][2]-(p[0][3][2]-p[0][1][2]))*one_div_dx2*one_div_dy3;
    b[39] = (p[3][3][2]-p[3][1][2]-(p[1][3][2]-p[1][1][2]))*one_div_dx3*one_div_dy3;

    // partial xz derivatives
    b[40] = (p[2][1][2]-p[2][1][0]-(p[0][1][2]-p[0][1][0]))*one_div_dx2*one_div_dz2;
    b[41] = (p[3][1][2]-p[3][1][0]-(p[1][1][2]-p[1][1][0]))*one_div_dx3*one_div_dz2;
    b[42] = (p[2][2][2]-p[2][2][0]-(p[0][2][2]-p[0][2][0]))*one_div_dx2*one_div_dz2;
    b[43] = (p[3][2][2]-p[3][2][0]-(p[1][2][2]-p[1][2][0]))*one_div_dx3*one_div_dz2;
    b[44] = (p[2][1][3]-p[2][1][1]-(p[0][1][3]-p[0][1][1]))*one_div_dx2*one_div_dz3;
    b[45] = (p[3][1][3]-p[3][1][1]-(p[1][1][3]-p[1][1][1]))*one_div_dx3*one_div_dz3;
    b[46] = (p[2][2][3]-p[2][2][1]-(p[0][2][3]-p[0][2][1]))*one_div_dx2*one_div_dz3;
    b[47] = (p[3][2][3]-p[3][2][1]-(p[1][2][3]-p[1][2][1]))*one_div_dx3*one_div_dz3;

    // partial yz derivatives
    b[48] = (p[1][2][2]-p[1][2][0]-(p[1][0][2]-p[1][0][0]))*one_div_dy2*one_div_dz2;
    b[49] = (p[2][2][2]-p[2][2][0]-(p[2][0][2]-p[2][0][0]))*one_div_dy2*one_div_dz2;
    b[50] = (p[1][3][2]-p[1][3][0]-(p[1][1][2]-p[1][1][0]))*one_div_dy3*one_div_dz2;
    b[51] = (p[2][3][2]-p[2][3][0]-(p[2][1][2]-p[2][1][0]))*one_div_dy3*one_div_dz2;
    b[52] = (p[1][2][3]-p[1][2][1]-(p[1][0][3]-p[1][0][1]))*one_div_dy2*one_div_dz3;
    b[53] = (p[2][2][3]-p[2][2][1]-(p[2][0][3]-p[2][0][1]))*one_div_dy2*one_div_dz3;
    b[54] = (p[1][3][3]-p[1][3][1]-(p[1][1][3]-p[1][1][1]))*one_div_dy3*one_div_dz3;
    b[55] = (p[2][3][3]-p[2][3][1]-(p[2][1][3]-p[2][1][1]))*one_div_dy3*one_div_dz3;

    // partial xyz derivatives
    b[56] = (p[2][2][2]-p[2][2][0]-(p[2][0][2]-p[2][0][0])-(p[0][2][2]-p[0][2][0])+(p[0][0][2]-p[0][0][0]))*one_div_dx2*one_div_dy2*one_div_dz2;
    b[57] = (p[3][2][2]-p[3][2][0]-(p[3][0][2]-p[3][0][0])-(p[1][2][2]-p[1][2][0])+(p[1][0][2]-p[1][0][0]))*one_div_dx3*one_div_dy2*one_div_dz2;
    b[58] = (p[2][3][2]-p[2][3][0]-(p[2][1][2]-p[2][1][0])-(p[0][3][2]-p[0][3][0])+(p[0][1][2]-p[0][1][0]))*one_div_dx2*one_div_dy3*one_div_dz2;
    b[59] = (p[3][3][2]-p[3][3][0]-(p[3][1][2]-p[3][1][0])-(p[1][3][2]-p[1][3][0])+(p[1][1][2]-p[1][1][0]))*one_div_dx3*one_div_dy3*one_div_dz2;
    b[60] = (p[2][2][3]-p[2][2][1]-(p[2][0][3]-p[2][0][1])-(p[0][2][3]-p[0][2][1])+(p[0][0][3]-p[0][0][1]))*one_div_dx2*one_div_dy2*one_div_dz3;
    b[61] = (p[3][2][3]-p[3][2][1]-(p[3][0][3]-p[3][0][1])-(p[1][2][3]-p[1][2][1])+(p[1][0][3]-p[1][0][1]))*one_div_dx3*one_div_dy2*one_div_dz3;
    b[62] = (p[2][3][3]-p[2][3][1]-(p[2][1][3]-p[2][1][1])-(p[0][3][3]-p[0][3][1])+(p[0][1][3]-p[0][1][1]))*one_div_dx2*one_div_dy3*one_div_dz3;
    b[63] = (p[3][3][3]-p[3][3][1]-(p[3][1][3]-p[3][1][1])-(p[1][3][3]-p[1][3][1])+(p[1][1][3]-p[1][1][1]))*one_div_dx3*one_div_dy3*one_div_dz3;

    for (uint i=0; i<64; i++)
    {
        double element = 0.0;
        for (uint j=0; j<64; j++)
            element += COEFF_MATRIX[i][j]*b[j];

        m_cached_coeffs[i_0*m_ny*m_nz+j_0*m_nz+k_0][i] = element;
    }
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
            m_x[i+1] = pos_of_grid_point(Dir::X, i, config); // +1 because we do not fill the two outermost elements yet

        for (uint j=0; j<m_ny; j++)
            m_y[j+1] = pos_of_grid_point(Dir::Y, j, config);

        for (uint k=0; k<m_nz; k++)
            m_z[k+1] = pos_of_grid_point(Dir::Z, k, config);

        set_grid_outermost();
    }
}


void Interpolator3D::set_grid_outermost()
{
    m_x[0] = m_x[1]-1.0;
    m_x[m_nx+1] = m_x[m_nx]+1.0;
    m_x[m_nx+2] = m_x[m_nx+1]+1.0;

    m_y[0] = m_y[1]-1.0;
    m_y[m_ny+1] = m_y[m_ny]+1.0;
    m_y[m_ny+2] = m_y[m_ny+1]+1.0;

    m_z[0] = m_z[1]-1.0;
    m_z[m_nz+1] = m_z[m_nz]+1.0;
    m_z[m_nz+2] = m_z[m_nz+1]+1.0;
}


void Interpolator3D::set_data_outermost()
{
    for (uint i=0; i<m_nx+3; i++)
        for (uint j=0; j<m_ny+3; j++)
            for (uint k=0; k<m_nz+3; k++)
                {
                    bool is_outermost = false;

                    uint i_temp = i;
                    if (i_temp==0)
                    {
                        i_temp = 1;
                        is_outermost = true;
                    }
                    else if (i_temp>m_nx)
                    {
                        i_temp = m_nx;
                        is_outermost = true;
                    }

                    uint j_temp = j;
                    if (j_temp==0)
                    {
                        j_temp = 1;
                        is_outermost = true;
                    }
                    else if (j_temp>m_ny)
                    {
                        j_temp = m_ny;
                        is_outermost = true;
                    }

                    uint k_temp = k;
                    if (k_temp==0)
                    {
                        k_temp = 1;
                        is_outermost = true;
                    }
                    else if (k_temp>m_nz)
                    {
                        k_temp = m_nz;
                        is_outermost = true;
                    }

                    if (is_outermost)
                        m_data[_INDEX(i, j, k)] = m_data[_INDEX(i_temp, j_temp, k_temp)];
                }
}


void Interpolator3D::export_data_old_format (const std::string& filepath) const
{
#ifdef _INTERP_LOG
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
#ifdef _INTERP_LOG
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
#ifdef _INTERP_LOG
    std::cout << "Finished exporting m_data to file" << std::endl;
#endif
}


void Interpolator3D::import_data_old_format (const std::string& filepath)
{
#ifdef _INTERP_LOG
    std::cout << "Importing m_data IN THE OLD FORMAT from file..." << std::endl;
#endif
    std::ifstream in;
    in.open(filepath);
    if (!in.is_open())
    {
#ifdef _INTERP_LOG
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
#ifdef _INTERP_LOG
        std::cerr << "Wrong m_data format, use generate_data() to generate m_data in the desired format. Aborting!" << std::endl;
#endif
        exit(54);
    }

    m_nx = std::stoi(line_vec[1]);
    m_ny = std::stoi(line_vec[2]);
    m_nz = std::stoi(line_vec[3]);
    
    set_grid(nullptr); // allocate memory for position arrays but leave the position arrays uninitialized

    prepare_data();

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
    set_data_outermost();

    cache_coeffs();
#ifdef _INTERP_LOG
    std::cout << "Finished importing m_data from file" << std::endl;
#endif
}


void Interpolator3D::export_data_plain_text (const std::string& filepath) const
{
#ifdef _INTERP_LOG
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
#ifdef _INTERP_LOG
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
#ifdef _INTERP_LOG
    std::cout << "Finished exporting m_data to file" << std::endl;
#endif
}


void Interpolator3D::import_data_plain_text (const std::string& filepath)
{
#ifdef _INTERP_LOG
    std::cout << "Importing m_data from file..." << std::endl;
#endif

    std::ifstream in;
    in.open(filepath);
    if (!in.is_open())
    {
#ifdef _INTERP_LOG
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
#ifdef _INTERP_LOG
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
#ifdef _INTERP_LOG
        std::cerr << "Wrong m_data format, use generate_data() to generate m_data in the desired format. Aborting!" << std::endl;
#endif
        exit(54);
    }

    m_nx = std::stoi(line_vec[1]);
    m_ny = std::stoi(line_vec[2]);
    m_nz = std::stoi(line_vec[3]);
    
    set_grid(nullptr); // allocate memory for position arrays but leave the position arrays uninitialized

    prepare_data();

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

    set_data_outermost();

    cache_coeffs();
#ifdef _INTERP_LOG
    std::cout << "Finished importing m_data from file" << std::endl;
#endif
}


void Interpolator3D::export_data (const std::string& filepath) const
{
    if (!m_data || !m_x || !m_y || !m_z)
        exit(57);

    std::ofstream out(filepath, std::ios::binary);
    if (!out.is_open())
        exit(52);

    // format check number
    out.write(reinterpret_cast<const char*>(&FORMAT_CHECK_NUM), sizeof(InterpFormatCheck));

    // nx, ny, nz
    out.write(reinterpret_cast<const char*>(&m_nx), sizeof(m_nx));
    out.write(reinterpret_cast<const char*>(&m_ny), sizeof(m_ny));
    out.write(reinterpret_cast<const char*>(&m_nz), sizeof(m_nz));

    // x, y, z
    out.write(reinterpret_cast<const char*>(m_x), sizeof(double)*(m_nx+3));
    out.write(reinterpret_cast<const char*>(m_y), sizeof(double)*(m_ny+3));
    out.write(reinterpret_cast<const char*>(m_z), sizeof(double)*(m_nz+3));

    // data
    out.write(reinterpret_cast<const char*>(m_data), sizeof(double)*(m_nx+3)*(m_ny+3)*(m_nz+3));

    out.close();
}


void Interpolator3D::import_data (const std::string& filepath)
{
    std::ifstream in(filepath, std::ios::binary);
    if (!in.is_open())
        exit(53);

    // check format of data to be imported
    InterpFormatCheck format_check_num;
    in.read(reinterpret_cast<char*>(&format_check_num), sizeof(InterpFormatCheck));

    if (format_check_num != FORMAT_CHECK_NUM)
        exit(54);

    // nx, ny, nz
    in.read(reinterpret_cast<char*>(&m_nx), sizeof(m_nx));
    in.read(reinterpret_cast<char*>(&m_ny), sizeof(m_ny));
    in.read(reinterpret_cast<char*>(&m_nz), sizeof(m_nz));

    // x, y, z
    set_grid(nullptr); // allocate pos arrays but leave them uninitialized
    in.read(reinterpret_cast<char*>(m_x), sizeof(double)*(m_nx+3));
    in.read(reinterpret_cast<char*>(m_y), sizeof(double)*(m_ny+3));
    in.read(reinterpret_cast<char*>(m_z), sizeof(double)*(m_nz+3));

    // data
    prepare_data();
    in.read(reinterpret_cast<char*>(m_data), sizeof(double)*(m_nx+3)*(m_ny+3)*(m_nz+3));

    in.close();

    cache_coeffs();
}


void Interpolator3D::generate_data (std::function<double (double,double,double)> func, const DataGenerationConfig* config, bool progress_monitor)
{
    set_grid(config);
    prepare_data();

    ProgressMonitor pm(m_nx);

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
#ifdef _INTERP_LOG
        if (progress_monitor)
        {
            pm.add_finished();
            pm.print_progress_percentage();
        }
#endif
    }

    set_data_outermost();

    cache_coeffs();
}


void Interpolator3D::find_closest_lower_data_point(int& i_0, int& j_0, int& k_0, double& x, double& y, double& z) const
{
    i_0 = find_index(m_x+1, m_nx, x);
    j_0 = find_index(m_y+1, m_ny, y);
    k_0 = find_index(m_z+1, m_nz, z);

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


double Interpolator3D::unicubic_interpolate (double p[4], double t[2], double z)
{
    double p2_m_p0_div_t2_m_t0 = (p[2]-p[0])/(1.0-t[0]);
    double p3_m_p1_div_t3_m_t1 = (p[3]-p[1])/t[1];
    double p1_m_p2 = p[1]-p[2];

    return p[1] + z*p2_m_p0_div_t2_m_t0 + z*z*(-3.0*p1_m_p2-2.0*p2_m_p0_div_t2_m_t0-p3_m_p1_div_t3_m_t1) + z*z*z*(2.0*p1_m_p2+p2_m_p0_div_t2_m_t0+p3_m_p1_div_t3_m_t1);
}


double Interpolator3D::bicubic_interpolate (double p[4][4], double t_y[2], double t_z[2], double y, double z)
{
	double unicubic_result[4];

	unicubic_result[0] = unicubic_interpolate(p[0], t_z, z);
	unicubic_result[1] = unicubic_interpolate(p[1], t_z, z);
	unicubic_result[2] = unicubic_interpolate(p[2], t_z, z);
	unicubic_result[3] = unicubic_interpolate(p[3], t_z, z);

	return unicubic_interpolate(unicubic_result, t_y, y);
}


double Interpolator3D::tricubic_interpolate (double p[4][4][4], double t_x[2], double t_y[2], double t_z[2], double x, double y, double z)
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
    t_x[1] = (t_x[3]-t_x[1])*one_div_t_x_2_t_x_1;

    t_y[0] = (t_y[0]-t_y[1])*one_div_t_y_2_t_y_1;
    t_y[1] = (t_y[3]-t_y[1])*one_div_t_y_2_t_y_1;

    double bicbuic_result_0 = bicubic_interpolate(p_z_0, t_x, t_y, x, y);
    double bicbuic_result_1 = bicubic_interpolate(p_z_1, t_x, t_y, x, y);

    return bicbuic_result_0+z*(bicbuic_result_1-bicbuic_result_0);
}


double Interpolator3D::get_interp_value_tricubic_old (double x, double y, double z) const
{
    int i_0, j_0, k_0;

    find_closest_lower_data_point(i_0, j_0, k_0, x, y, z);

    double p[4][4][4];

    for (int i=0; i<4; i++)
        for (int j=0; j<4; j++)
            for (int k=0; k<4; k++)
                p[i][j][k] = m_data[_INDEX(i+i_0, j+j_0, k+k_0)];

    double t_x[2];
    double t[4];
    for (int i=0; i<4; i++)
        t[i] = m_x[i+i_0];
    double one_div_t_2_m_t_1 = 1.0/(t[2]-t[1]);
    x = (x-t[1])*one_div_t_2_m_t_1;
    t_x[0] = (t[0]-t[1])*one_div_t_2_m_t_1;
    t_x[1] = (t[3]-t[1])*one_div_t_2_m_t_1;

    double t_y[2];
    for (int i=0; i<4; i++)
        t[i] = m_y[i+j_0];
    double one_div_t_2_m_t_1 = 1.0/(t[2]-t[1]);
    y = (y-t[1])*one_div_t_2_m_t_1;
    t_y[0] = (t[0]-t[1])*one_div_t_2_m_t_1;
    t_y[1] = (t[3]-t[1])*one_div_t_2_m_t_1;

    double t_z[2];
    for (int i=0; i<4; i++)
        t[i] = m_z[i+k_0];
    double one_div_t_2_m_t_1 = 1.0/(t[2]-t[1]);
    z = (z-t[1])*one_div_t_2_m_t_1;
    t_z[0] = (t[0]-t[1])*one_div_t_2_m_t_1;
    t_z[1] = (t[3]-t[1])*one_div_t_2_m_t_1;

    return tricubic_interpolate(p, t_x, t_y, t_z, x, y, z);
}


inline double Interpolator3D::A(Coeffs& coeffs, uint i, uint j, double z, double z2, double z3)
{
    return coeffs[i+j*4+0*16] + coeffs[i+j*4+1*16]*z + coeffs[i+j*4+2*16]*z2 + coeffs[i+j*4+3*16]*z3;
}


inline double Interpolator3D::B(Coeffs& coeffs, uint i, double y, double y2, double y3, double z, double z2, double z3)
{
    return A(coeffs, i, 0, z, z2, z3) + A(coeffs, i, 1, z, z2, z3)*y + A(coeffs, i, 2, z, z2, z3)*y2 + A(coeffs, i, 3, z, z2, z3)*y3;
}


double Interpolator3D::get_interp_value_tricubic(double x, double y, double z) const
{
    int i_0, j_0, k_0;

    find_closest_lower_data_point(i_0, j_0, k_0, x, y, z);

    x = (x-m_x[i_0+1])/(m_x[i_0+2]-m_x[i_0+1]);
    y = (y-m_y[j_0+1])/(m_y[j_0+2]-m_y[j_0+1]);
    z = (z-m_z[k_0+1])/(m_z[k_0+2]-m_z[k_0+1]);

    double x2 = x*x;
    double x3 = x2*x;

    double y2 = y*y;
    double y3 = y2*y;

    double z2 = z*z;
    double z3 = z2*z;

    return B(m_cached_coeffs[i_0*m_ny*m_nz + j_0*m_nz + k_0], 0, y, y2, y3, z, z2, z3)
        + B(m_cached_coeffs[i_0*m_ny*m_nz + j_0*m_nz + k_0], 1, y, y2, y3, z, z2, z3)*x
        + B(m_cached_coeffs[i_0*m_ny*m_nz + j_0*m_nz + k_0], 2, y, y2, y3, z, z2, z3)*x2
        + B(m_cached_coeffs[i_0*m_ny*m_nz + j_0*m_nz + k_0], 3, y, y2, y3, z, z2, z3)*x3;
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
    prepare_data();
    set_grid(nullptr);

    prepare_cached_coeffs();

    std::copy(other.m_data, other.m_data+(m_nx+3)*(m_ny+3)*(m_nz+3), m_data);
    std::copy(other.m_cached_coeffs, other.m_cached_coeffs+m_nx*m_ny*m_nz, m_cached_coeffs);
    std::copy(other.m_x, other.m_x+m_nx+3, m_x);
    std::copy(other.m_y, other.m_y+m_ny+3, m_y);
    std::copy(other.m_z, other.m_z+m_nz+3, m_z);

}


Interpolator3D::Interpolator3D (Interpolator3D&& other)
    : m_nx(other.m_nx)
    , m_ny(other.m_ny)
    , m_nz(other.m_nz)
    , m_data(other.m_data)
    , m_cached_coeffs(other.m_cached_coeffs)
    , m_x(other.m_x)
    , m_y(other.m_y)
    , m_z(other.m_z)
{
    other.m_data = nullptr;
    other.m_cached_coeffs = nullptr;
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

    prepare_data();
    prepare_cached_coeffs();
    set_grid(nullptr);

    std::copy(other.m_data, other.m_data+(m_nx+3)*(m_ny+3)*(m_nz+3), m_data);
    std::copy(other.m_cached_coeffs, other.m_cached_coeffs+m_nx*m_ny*m_nz, m_cached_coeffs);
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
    other.m_cached_coeffs = nullptr;
    other.m_x = nullptr;
    other.m_y = nullptr;
    other.m_z = nullptr;

    return *this;
}


Interpolator3D::~Interpolator3D()
{
    safe_delete_grid();
    
    safe_delete_data();

    safe_delete_cached_coeffs();
}


#undef _INDEX