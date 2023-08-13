#pragma once


#include <stdlib.h>
#include <string>


enum class Dir
{
    x, y, z
};


struct IndicesVec
{
    int i, j, k;
};


struct DataGenerationConfig
{
    int n_x; double x_min; double x_max; std::string x_grid_spacing = "linear";
    int n_y; double y_min; double y_max; std::string y_grid_spacing = "linear";
    int n_z; double z_min; double z_max; std::string z_grid_spacing = "linear";
};


class Interpolator3D
{
public:

    double* data_array;

    bool data_array_is_prepared = false;

    double* x_pos;
    double* y_pos;
    double* z_pos;

    bool grid_is_set = false;

    unsigned int n_x = 0;
    unsigned int n_y = 0;
    unsigned int n_z = 0;

    double safe_get_x_pos(int i);

    double safe_get_y_pos(int i);

    double safe_get_z_pos(int i);

private:

    void set_grid(DataGenerationConfig* config);

    void delete_grid();

    void delete_data_array();

    void prepare_data_array();

    uint index(uint i, uint j, uint k);    

    double pos_of_grid_point(Dir dir, int index, DataGenerationConfig* config);

    IndicesVec find_indices_of_closest_lower_data_point(double x, double y, double z);

public:

    void generate_data(double func(double x, double y, double z), DataGenerationConfig* config, bool progress_monitor);

    void export_data(std::string filepath);

    void import_data(std::string filepath);

    double safe_get_data_point(int i, int j, int k);

    double unicubic_interpolate(double p[4], double t_z[4], double z);

    double bicubic_interpolate(double p[4][4], double t_y[4], double t_z[4], double y, double z);

    double tricubic_interpolate(double p[4][4][4], double t_x[4], double t_y[4], double t_z[4], double x, double y, double z);

    double get_interp_value_tricubic(double x, double y, double z);

    void setup_gsl_interp();

    double get_interp_value_gsl_tricubic(double x, double y, double z);

    ~Interpolator3D();
};