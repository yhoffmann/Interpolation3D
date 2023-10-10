#pragma once


#include <stdlib.h>
#include <string>
#include <math.h>
#include <functional>
#include <mutex>
#include <omp.h>
#include <iostream>


enum class Dir : unsigned char
{
    x, y, z
};


enum GridSpacing : unsigned char
{
    Linear, Exponential
};


struct DataGenerationConfig
{
    uint n_x=300; double x_min=0.0; double x_max=15.0; GridSpacing x_grid_spacing = Exponential; double x_exp_grid_spacing_parameter = 8.0;
    uint n_y=300; double y_min=0.0; double y_max=15.0; GridSpacing y_grid_spacing = Exponential; double y_exp_grid_spacing_parameter = 8.0;
    uint n_z=30; double z_min=0.0; double z_max=M_PI; GridSpacing z_grid_spacing = Linear; double z_exp_grid_spacing_parameter = 8.0;
};


class Interpolator3D
{public:
    double* data_array = nullptr;
    double* x_pos = nullptr;
    double* y_pos = nullptr;
    double* z_pos = nullptr;
    uint n_x = 0;
    uint n_y = 0;
    uint n_z = 0;

    void safe_delete_grid();
    void set_grid(const DataGenerationConfig* config);
    void set_grid_outermost();
    void set_data_array_outermost();
    void safe_delete_data_array();
    void prepare_data_array();
    double pos_of_grid_point(Dir dir, int index, const DataGenerationConfig* config) const;
    void find_closest_lower_data_point(int& i_0, int& j_0, int& k_0, double& x, double& y, double& z) const;
    static uint find_index(double* arr, int size, double x);

public:

    void generate_data(std::function<double (double,double,double)> func, const DataGenerationConfig* config, bool progress_monitor);
    void export_data_old_format(const std::string& filepath) const;
    void import_data_old_format(const std::string& filepath);
    void export_data (const std::string& filepath) const;
    void import_data(const std::string& filepath);
    double safe_get_data_point(int i, int j, int k) const;
    double safe_get_x_pos(int i) const;
    double safe_get_y_pos(int i) const;
    double safe_get_z_pos(int i) const;
    static double unicubic_interpolate(double p[4], double t_z[4], double z);
    static double bicubic_interpolate(double p[4][4], double t_y[4], double t_z[4], double y, double z);
    static double tricubic_interpolate(double p[4][4][4], double t_x[4], double t_y[4], double t_z[4], double x, double y, double z);
    double get_interp_value_tricubic(double x, double y, double z) const;
    double get_interp_value_bicubic_unilinear(double x, double y, double z) const;
    
    Interpolator3D();
    Interpolator3D(const std::string& filepath);
    Interpolator3D(const Interpolator3D& other);
    Interpolator3D& operator=(const Interpolator3D& other);
    Interpolator3D& operator=(Interpolator3D&& other);
    ~Interpolator3D();
};

inline std::mutex find_index_lock;