#include "../include/TricubicInterpolator.h"
#include <cmath>


XYZ TricubicInterpolator::inputs_for_grid_pos(luint i, luint j, luint k, DataGenerationConfig& config)
{
    XYZ xyz;

    // log only really makes sense for values > 0
    if (config.x_grid_spacing == "linear")
    {
        xyz.x = config.x_min+(config.x_max-config.x_min)*(double(i))/(double(config.n_x-1));
    }
    else if (config.x_grid_spacing == "log")
    {
        if (config.x_min == 0) { config.x_min = 1.0e-3*(config.x_max-config.x_min)/(double(config.n_x-1)); }
        xyz.x = config.x_min*exp(log(config.x_max/config.x_min)*(double(i))/double(config.n_x-1))-config.n_x;
    }

    if (config.x_grid_spacing == "linear")
    {
        xyz.y = config.y_min+(config.y_max-config.y_min)*(double(i))/(double(config.n_y));
    }
    else if (config.y_grid_spacing == "log")
    {
        if (config.y_min == 0) { config.y_min = 1.0e-3*(config.y_max-config.y_min)/(double(config.n_y-1)); }
        xyz.y = config.y_min*exp(log(config.y_max/config.y_min)*(double(i))/double(config.n_y-1))-config.n_y;
    }

    if (config.x_grid_spacing == "linear")
    {
        xyz.z = config.z_min+(config.z_max-config.z_min)*(double(i))/(double(config.n_z));
    }
    else if (config.z_grid_spacing == "log")
    {
        if (config.z_min == 0) { config.z_min = 1.0e-3*(config.z_max-config.z_min)/(double(config.n_z-1)); }
        xyz.z = config.z_min*exp(log(config.z_max/config.z_min)*(double(i))/double(config.n_z-1))-config.n_z;
    }

    return xyz;
}


void TricubicInterpolator::set_grid_pos (DataGenerationConfig config)
{
    luint n_x = config.n_x;
    luint n_y = config.n_y;
    luint n_z = config.n_z;
    
    // initializing 3d vector and resizing it to size n_x,n_y,n_z
    grid_pos.resize(n_x);
    for (luint i=0; i<n_x; i++)
    {
        grid_pos[i].resize(n_y);
    }
    for (luint i=0; i<n_x; i++)
    {
        for (luint j=0; j<n_y; j++)
        {
            grid_pos[i][j].resize(n_z);
        }
    }

    for (luint i=0; i<n_x; i++)
    {
        for (luint j=0; j<n_y; j++)
        {
            for (luint k=0; k<n_z; k++)
            {
                grid_pos[i][j][k] = inputs_for_grid_pos(i,j,k,config);
            }
        }
    }
}


void TricubicInterpolator::generate_data (double func(double x, double y, double z), DataGenerationConfig config, bool progress_monitor)
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
                new_data[i][j][k] = func(grid_pos[i][j][k].x,grid_pos[i][j][k].y,grid_pos[i][j][k].z);
            }
        }
    }
}