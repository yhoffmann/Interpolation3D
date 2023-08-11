#include "../include/Interpolator3D.h"
#include <cmath>
#include <iostream>
#include <iomanip>


double f(double x, double y, double z)
{
    return log(x+1.0)*log(y+2.0)*exp(-z*z);
}


int main (int argc, char** argv)
{
    std::string filepath = "InterpolatorData/datatest.txt";

    Interpolator3D ip;

    DataGenerationConfig conf;

    conf.n_x = 101;
    conf.n_y = 101;
    conf.n_z = 101;

    //conf.x_grid_spacing = "log";
    //conf.y_grid_spacing = "log";
    //conf.z_grid_spacing = "log";

    conf.x_min = 1.0;
    conf.y_min = 1.0;
    conf.z_min = 1.0;

    conf.x_max = 11;
    conf.y_max = 11;
    conf.z_max = 11;

    ip.generate_data(f,conf,false);
    //ip.load_data(filepath);
/*
    for (int i=0; i<200; i++)
    {
        double x = -10.0+20.0*double(i)/double(200)-0.01;
        (void)ip.tricubic_get_value(x,0,0);
        //std::cout << ip.safe_get_pos(X,i) << " " << ip.safe_get_data_point(i,0,0) << std::endl;
    }
*/

    double y = 1.0;
    double z = 1.0;
    double x = 1.0;
    //std::cout << x << " " << f(x,y,z) << " " << ip.tricubic_get_value_nonreg(x,y,z) << std::endl;
    int imax = 1e8+1;
    for (int i=0; i<imax; i++)
    {
        double x = 5.0*double(i)/double(imax-1);
        //std::cout << std::setprecision(10) << x << " " << f(x,y,z) << " " << ip.tricubic_get_value(x,y,z) << std::endl;
        (void)ip.tricubic_get_value(x,y,z);
    }

    //ip.print_data_to_file(filepath);

    return 0;
}