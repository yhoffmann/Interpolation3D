#include "../include/Interpolator3D.hpp"
// #include "../../utils/cubature/cubature.h"
#include <chrono>

double func(double x, double y, double z) {
  return log(x * x * y * y + 1.2) * (1.0 + 0.1 * cos(2.0 * y)) +
         exp(-z * (x - 0.1) / (y + 0.1)) / (1.0 + 1.0 / (x + 1.0));
}

// Interpolator a;
// int inte (unsigned ndim, const double* xx, void* userdata, unsigned fdim,
// double* ff)
// {

//     ff[0] = a(xx[0], xx[1], xx[2]);
//     *(uint*)userdata += 1;
//     return 0;
// }

// void integrate()
// {
//     double xmin[3];
//     double xmax[3];
//     uint count = 0;
//     xmin[0] = 0.0; xmax[0] = 15.0;
//     xmin[1] = 0.0; xmax[1] = 15.0;
//     xmin[2] = 0.0; xmax[2] = 15.0;

//     double ret(0), err(0);

//     hcubature(1, inte, &count, 3, xmin, xmax, 1e8, 0.0, 0.0,
//     error_norm::ERROR_INDIVIDUAL, &ret, &err);

//     std::cout << ret << " " << err << " " << count << std::endl;
// }

int main() {
  Interpolator3D inter;
  inter.set_num_threads(32);

  Interpolator3D::DataGenerationConfig c;
  // c.nx = 300+1;
  // c.ny = 300+1;
  // c.x_exp_grid_spacing_parameter = 6.0;
  // c.y_exp_grid_spacing_parameter = 6.0;
  c.z_max = 15.0;

  inter.generate_data(func, &c, false);

  std::cout << "starting" << std::endl;

  auto start = std::chrono::high_resolution_clock::now();

  double sum = 0.0;
  uint num_runs = 3e8;
  for (uint i = 0, imax = num_runs; i < imax; ++i) {
    double x = double(i) / double(imax - 1) * 15.0;
    double y = double(i) / double(imax - 1) * 8.0; // 3*15.0/50.0;
    double z = -double(i) / double(imax - 1) * 4.0 + 7.0; //*double(i)/double(imax-1)*15.0;//double(i)/double(imax-1)*4.0
                                                          //+ 1.0e-6;

    volatile double val = inter(x, y, z, Interpolator3D::Tricubic);
    // sum += val;
    // volatile double val = func(x, y, z);

    // volatile double val_old = inter.get_interp_value_tricubic_old(x, y, z);
    // volatile double val_file = inter2.get_interp_value_tricubic_old(x, y, z);

    // std::cout x << " " << func(x, y, z) << " " << val << " " <<
    // std::abs(func(x, y, z)-val)/func(x, y, z) << std::endl; std::cout << x <<
    // " " << val << " " << val_old << " " << std::abs(val_old-val)/val <<
    // std::endl; std::cout << x << " " << val << " " << val_file << " " <<
    // std::abs(val_file-val)/val << std::endl;

    // for (uint I=0; I<4; I++)
    //     for (uint J=0; J<4; J++)
    //         for (uint K=0; K<4; K++)
    //             sum += coeffs[I*16+J*4+K] * pow(x, I) * pow(y, J) * pow(z,
    //             K);

    // sum -= 1.0;
  }
  // std::cout << "Sum is " << sum << std::endl;

  auto end = std::chrono::high_resolution_clock::now();
  double time_ms =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
          .count();

  std::cout << "Just the interpolation part loop took " << time_ms << "ms"
            << std::endl;
  std::cout << "That results in " << double(num_runs) / time_ms * 1.0e3 * 1.0e-6
            << " million interpolations per second" << std::endl;

  return 0;
}
