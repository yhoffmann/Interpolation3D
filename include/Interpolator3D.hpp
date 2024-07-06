#pragma once


#include <stdlib.h>
#include <string>
#include <math.h>
#include <functional>


enum class Dir : unsigned char
{
    X, Y, Z
};


enum GridSpacing : unsigned char
{
    Linear, Exponential
};


struct DataGenerationConfig
{
    uint nx=300; double x_min=0.0; double x_max=15.0; GridSpacing x_grid_spacing = Exponential; double x_exp_grid_spacing_parameter = 8.0;
    uint ny=300; double y_min=0.0; double y_max=15.0; GridSpacing y_grid_spacing = Exponential; double y_exp_grid_spacing_parameter = 8.0;
    uint nz=30; double z_min=0.0; double z_max=M_PI; GridSpacing z_grid_spacing = Linear; double z_exp_grid_spacing_parameter = 8.0;
};


enum InterpolationType
{
    BicubicUnilinear, Tricubic
};


struct Coeffs
{
private:
    double m_coeffs[64];

public:
    inline double& operator[](uint i)
    {
        return m_coeffs[i];
    }
};


class Interpolator3D
{
protected:

    uint m_nx = 0;
    uint m_ny = 0;
    uint m_nz = 0;
    double* m_data = nullptr;
    Coeffs* m_cached_coeffs = nullptr;
    double* m_x = nullptr;
    double* m_y = nullptr;
    double* m_z = nullptr;

public:

    static double step_growth_factor(double k);
    void generate_data(std::function<double (double,double,double)> func, const DataGenerationConfig* config, bool progress_monitor);
    void export_data_old_format(const std::string& filepath) const;
    void import_data_old_format(const std::string& filepath);
    void export_data_plain_text(const std::string& filepath) const;
    void import_data_plain_text(const std::string& filepath);
    void export_data(const std::string& filepath) const;
    void import_data(const std::string& filepath);
    static double unicubic_interpolate(double p[4], double t_z[2], double z);
    static double bicubic_interpolate(double p[4][4], double t_y[2], double t_z[2], double y, double z);
    static double tricubic_interpolate(double p[4][4][4], double t_x[2], double t_y[2], double t_z[2], double x, double y, double z);
    double get_interp_value_tricubic(double x, double y, double z) const;
    double get_interp_value_tricubic_old(double x, double y, double z) const;
    double get_interp_value_bicubic_unilinear(double x, double y, double z) const;
    
    constexpr double operator() (double x, double y, double z, InterpolationType type = BicubicUnilinear) const
    {
        switch (type)
        {
            case BicubicUnilinear:
                return get_interp_value_bicubic_unilinear(x, y, z);
            break;

            case Tricubic: // fallthrough

            default:
                return get_interp_value_tricubic(x, y, z);
            break;
        }
    }

    Interpolator3D();
    Interpolator3D(const std::string& filepath);
    Interpolator3D(const Interpolator3D&);
    Interpolator3D(Interpolator3D&&);
    Interpolator3D& operator=(const Interpolator3D&);
    Interpolator3D& operator=(Interpolator3D&&);
    ~Interpolator3D();

protected:

    void safe_delete_grid();
    void set_grid(const DataGenerationConfig* config);
    void set_grid_outermost();
    void set_data_outermost();
    void safe_delete_data();
    void prepare_data();
    void prepare_cached_coeffs();
    void safe_delete_cached_coeffs();
    void cache_coeffs();
    void set_single_cell_coeffs(uint i_0, uint j_0, uint k_0);
    static double A(Coeffs& coeffs, uint i, uint j, double z, double z2, double z3);
    static double B(Coeffs& coeffs, uint i, double y, double y2, double y3, double z, double z2, double z3);
    double pos_of_grid_point(Dir dir, int index, const DataGenerationConfig* config) const;
    void find_closest_lower_data_point(int& i_0, int& j_0, int& k_0, double& x, double& y, double& z) const;
    static uint find_index(double* arr, int size, double x);
};