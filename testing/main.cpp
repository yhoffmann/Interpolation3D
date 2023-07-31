#include "../include/Interpolator3D.h"
#include <cmath>
#include <iostream>
#include <iomanip>


double f(double x, double y, double z)
{
    return x*x*y*y*cos(z)*exp(-x*x-y*y-2.0*x*y*cos(z));
}


int main (int argc, char** argv)
{
    std::string filepath = "InterpolatorData/datatest.txt";

    Interpolator3D ip;

    DataGenerationConfig conf;

    conf.n_x = 1000;
    conf.n_y = 1000;
    conf.n_z = 1000;

    //conf.x_grid_spacing = "log";
    //conf.y_grid_spacing = "log";
    //conf.z_grid_spacing = "log";

    conf.x_min = 0;
    conf.y_min = 0;
    conf.z_min = -0.1;

    conf.x_max = 5;
    conf.y_max = 5;
    conf.z_max = M_PI+0.1;

    ip.generate_data(f,conf,false);
    //ip.load_data(filepath);
/*
    double p[4] = {-0.4, 1.0, 2.5, 1.5};
    for (int i=0; i<100; i++)
    {
        double x = -1.0+3.0*i/99.0;
        std::cout << x << " " << ip.unicubic_interpolate(p,x) << std::endl;
    }
*/
    double x = 1.5;
    double y = 1.0;
    double z = 1.0;
    std::cout << std::setprecision(10) << ip.x_pos[2] << " " << ip.y_pos[1] << " " << ip.z_pos[0] << " " << ip.data_array[2][1][0] << " " << f(1.0,0.5,-0.1) << std::endl;
    std::cout << ip.trilinear_get_value(x,y,z) << " " << f(x,y,z) << ", err: " << (ip.trilinear_get_value(x,y,z)-f(x,y,z))/f(x,y,z) << std::endl;

    std::cout << ip.bicubic_unilinear_get_value(x,y,z) << " " << f(x,y,z) << ", err: " << (ip.bicubic_unilinear_get_value(x,y,z)-f(x,y,z))/f(x,y,z) << std::endl;

    std::cout << ip.tricubic_get_value(x,y,z) << " " << f(x,y,z) << ", err: " << (ip.tricubic_get_value(x,y,z)-f(x,y,z))/f(x,y,z) << std::endl;
/*
    luint imax = 551;
    for (luint i=0; i<imax; i++)
    {
        double x = -0.5+(5-(-0.5))*double(i)/double(imax-1);
        std::cout << x << " " << ip.trilinear_get_value(x,y,z) << " " << ip.bicubic_unilinear_get_value(x,y,z) << " " << ip.tricubic_get_value(x,y,z) << " " << f(x,y,z) << " " << ip.bicubic_unilinear_get_value(x,y,z) << std::endl;
        //(void)ip.tricubic_get_value_test(x,y,z);
        //(void)ip.bicubic_unilinear_get_value(x,y,z);
    }
*/
    //ip.print_data_to_file(filepath);
    
    return 0;
}