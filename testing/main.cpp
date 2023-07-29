#include "../include/Interpolator3D.h"
#include <cmath>
#include <iostream>


double f(double x, double y, double z)
{
    return exp(-x*x-y*y-z*z);
}


int main (int argc, char** argv)
{
    std::string filepath = "InterpolatorData/datatest.txt";

    Interpolator3D ip;

    DataGenerationConfig conf;

    conf.n_x = 100;
    conf.n_y = 100;
    conf.n_z = 100;

    conf.x_grid_spacing = "linear";
    conf.y_grid_spacing = "linear";
    conf.z_grid_spacing = "linear";

    conf.x_min = 0;
    conf.y_min = 0;
    conf.z_min = 0;

    conf.x_max = 5;
    conf.y_max = 5;
    conf.z_max = 5;

    ip.generate_data(f,conf,false);
    //ip.load_data(filepath);

    double x = 3;
    double y = 0.5;
    double z = 0.5;
    
    std::cout << ip.trilinear_get_value(x,y,z) << " " << f(x,y,z) << ", err: " << (ip.trilinear_get_value(x,y,z)-f(x,y,z))/f(x,y,z) << std::endl;

    std::cout << ip.bicubic_unilinear_get_value(x,y,z) << " " << f(x,y,z) << ", err: " << (ip.bicubic_unilinear_get_value(x,y,z)-f(x,y,z))/f(x,y,z) << std::endl;
    

    //ip.print_data_to_file(filepath);
    
    return 0;
}