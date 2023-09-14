#include "../include/Interpolator3D.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <thread>


double f(double x, double y, double z)
{
    return log(x+1.0)*log(y+2.0)*exp(-z*z);
}


int main (int argc, char** argv)
{
    std::string filepath = "InterpolatorData/datatest.txt";

    Interpolator3D ip;

    DataGenerationConfig conf;

    conf.n_x = 200;
    conf.n_y = 200;
    conf.n_z = 140;

    //conf.x_grid_spacing = "log";
    //conf.y_grid_spacing = "log";
    //conf.z_grid_spacing = "log";

    conf.x_min = 0.0;
    conf.y_min = 0.0;
    conf.z_min = 0.0;

    conf.x_max = 10;
    conf.y_max = 10;
    conf.z_max = 10;

    ip.generate_data(f,&conf,false);
    //ip.load_data(filepath);
/*
    for (int i=0; i<200; i++)
    {
        double x = -10.0+20.0*double(i)/double(200)-0.01;
        (void)ip.tricubic_get_value(x,0,0);
        //std::cout << ip.safe_get_pos(X,i) << " " << ip.safe_get_data_point(i,0,0) << std::endl;
    }
*/

    double x = 5.3;
    double y = 2.158;
    double z = 1.012;
    std::cout << x << " " << f(x,y,z) << "\n" << ip.get_interp_value_tricubic(x,y,z) << std::endl;

    int imax = 1e8+1;

    double sum = 0.0;
    double current_value = 0.0;
    //#pragma omp parallel for
    for (int i=0; i<imax; i++)
    {
        double x = 10.0*double(i)/double(imax-1)+0.01;
        double y=x, z=x;
        //std::cout << std::setprecision(10) << x << " " << f(x,y,z) << " " << (f(x,y,z)-ip.get_interp_value_tricubic(x,y,z))/f(x,y,z) << std::endl;
        (void)ip.get_interp_value_tricubic(x,y,z);
        //(void)ip.get_interp_value_bicubic_unilinear(x,y,z);
        //current_value = (ip.get_interp_value_bicubic_unilinear(x,y,z)-f(x,y,z))/f(x,y,z);
        //sum += current_value*current_value;
        //std::cout << "x " << x << std::endl;
    }
    //std::cout << std::sqrt(sum) << std::endl;

    //ip.print_data_to_file(filepath);
}