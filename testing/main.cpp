#include "../include/Interpolator3D.h"
#include <cmath>
#include <iostream>
#include <iomanip>


double f(double x, double y, double z)
{
    return log(x+1.0)*exp(-x*y*z-y*y/2.0-cos(z)*cos(z));
}


void g(double x[2])
{
    std::cout << ++x[1] << std::endl;
}


int main (int argc, char** argv)
{
    std::string filepath = "InterpolatorData/datatest.txt";

    Interpolator3D ip;

    DataGenerationConfig conf;

    conf.n_x = 201;
    conf.n_y = 201;
    conf.n_z = 141;

    conf.x_grid_spacing = "log";
    conf.y_grid_spacing = "log";
    //conf.z_grid_spacing = "log";

    conf.x_min = 0.0;
    conf.y_min = 0.0;
    conf.z_min = 0.0;

    conf.x_max = 10;
    conf.y_max = 10;
    conf.z_max = 10;

    ip.generate_data(f,conf,false);
    //ip.load_data(filepath);

    double y = 1.0;
    double z = 1.0;
    double x = 1.0;
    //std::cout << x << " " << f(x,y,z) << " " << ip.tricubic_get_value_nonreg(x,y,z) << std::endl;
    int imax = 1e3+1;
    for (int i=0; i<imax; i++)
    {
        double x = 1.0*double(i)/double(imax-1);
        std::cout << std::setprecision(10) << x << " " << f(x,y,z) << " " << ip.tricubic_get_value_nonreg(x,y,z) << std::endl;
    }

    //ip.print_data_to_file(filepath);

    return 0;
}