#pragma once


#include <stdlib.h>
#include <vector>
#include <string>


typedef std::vector<std::vector<std::vector<double>>> vec_3d;


enum XYZ
{
    X, Y, Z
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

    vec_3d data_array;
    std::vector<double> x_pos;
    std::vector<double> y_pos;
    std::vector<double> z_pos;

    unsigned int n_x;
    unsigned int n_y;
    unsigned int n_z;

    double safe_get_x_pos(int i);

    double safe_get_y_pos(int i);

    double safe_get_z_pos(int i);

    void set_grid(DataGenerationConfig& config);

    void prepare_data_array();

    double pos_of_grid_point(XYZ xyz, int index, DataGenerationConfig& config);

    IndicesVec find_indices_of_closest_lower_data_point(double x, double y, double z);


public:

    void generate_data(double func(double x, double y, double z), DataGenerationConfig& config, bool progress_monitor);

    void export_data(std::string filepath);

    void import_data(std::string filepath);

    double safe_get_data_point(int i, int j, int k);

    double unicubic_interpolate(double p[4], double t_z[4], double z);

    double bicubic_interpolate(double p[4][4], double t_y[4], double t_z[4], double y, double z);

    double tricubic_interpolate(double p[4][4][4], double t_x[4], double t_y[4], double t_z[4], double x, double y, double z);

    double tricubic_get_value(double x, double y, double z);
};