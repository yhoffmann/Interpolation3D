#pragma once


#include <stdlib.h>
#include <vector>
#include <string>
#include <iostream>

typedef long unsigned int luint;


struct XYZ
{
    double x, y, z;
};


typedef std::vector<std::vector<std::vector<double>>> vec_3d;
typedef std::vector<std::vector<std::vector<XYZ>>> vec_grid_positions;


struct DataGenerationConfig
{
    luint n_x; double x_min; double x_max; std::string x_grid_spacing;
    luint n_y; double y_min; double y_max; std::string y_grid_spacing;
    luint n_z; double z_min; double z_max; std::string z_grid_spacing;
};


class TricubicInterpolator
{
private:

    vec_3d data;
    vec_grid_positions grid_pos;

    bool data_is_loaded = false;
    bool grid_pos_is_set = false;

    void set_grid_pos(DataGenerationConfig config);

    bool print_data_to_file();

    XYZ inputs_for_grid_pos(luint i, luint j, luint k, DataGenerationConfig& config);


public:

    TricubicInterpolator();

    void generate_data(double func(double x, double y, double z), DataGenerationConfig config, bool progress_monitor);

    bool load_data(double func(double x, double y, double z));

    double get_value(double x, double y, double z);
};