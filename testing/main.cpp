#include "../include/Interpolator3D.h"
#include <cmath>
#include <iostream>
#include <iomanip>


double f(double x, double y, double z)
{
    return log(x+1.0)*exp(-y*y-z*z);
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

    conf.n_x = 6;
    conf.n_y = 4;
    conf.n_z = 4;

    conf.x_grid_spacing = "log";
    //conf.y_grid_spacing = "log";
    //conf.z_grid_spacing = "linear";

    conf.x_min = 0.0;
    conf.y_min = 0.0;
    conf.z_min = 0.0;

    conf.x_max = 10;
    conf.y_max = 4;
    conf.z_max = 4;

    ip.generate_data(f,conf,false);
    //ip.load_data(filepath);

    ip.setup_gsl_interp();

    double t[4] = {0.0, 1.0, 2.0, 3.0};
    double p[4] = {0.0, 1.5, 0.2, 0.5};

    for (int i=0; i<100; i++)
    {
        double x = 1.0*i/99.0;
        std::cout << x << " " << ip.unicubic_interpolate_nonreg(p,t,x) << " " << ip.unicubic_interpolate(p,x) << std::endl;
    }

    

    //ip.print_data_to_file(filepath);

    
    
    return 0;
}