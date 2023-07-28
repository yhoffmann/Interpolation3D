#pragma once


#include <stdlib.h>
#include <vector>
#include <string>
#include <iostream>

typedef long unsigned int luint;


enum XYZ
{
    x, y, z
};


typedef std::vector<std::vector<std::vector<double>>> vec_3d;


struct DataGenerationConfig
{
    luint n_x; double x_min; double x_max; std::string x_grid_spacing = "linear";
    luint n_y; double y_min; double y_max; std::string y_grid_spacing = "linear";
    luint n_z; double z_min; double z_max; std::string z_grid_spacing = "linear";
};


class Interpolator3D
{
public:

    vec_3d data;
    std::vector<double> x_pos;
    std::vector<double> y_pos;
    std::vector<double> z_pos;

    bool data_is_loaded = false;
    bool grid_pos_is_set = false;

    void set_grid_pos(DataGenerationConfig& config);

    void print_data_to_file(std::string filepath);

    double inputs_for_pos(XYZ xyz, luint index, DataGenerationConfig& config);


public:

    //TrilinearInterpolator();

    void generate_data(double func(double x, double y, double z), DataGenerationConfig& config, bool progress_monitor);

    void load_data(std::string filepath);

    double trilinear_get_value(double x, double y, double z);

    double tricubic_get_value(double x, double y, double z);
};