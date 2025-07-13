# Interpolation3D

This class allows for tricubic interpolation with arbitrary grid spacing.

The tricubic interpolation is based on [this paper](http://www.cds.caltech.edu/~marsden/bib/2005/08-LeMa2005/LeMa2005.pdf)
>Lekien, F. and Marsden, J. (2005), Tricubic interpolation in three dimensions. Int. J. Numer. Meth. Engng., 63: 455-471. https://doi.org/10.1002/nme.1296  

by Lekien and Marsden. Unfortunately, the original source code by the authors does not seem be available anymore (as also described [here](https://github.com/nbigaouette/libtricubic/), which is an unofficial mirror of version `1.0` of the original authors' implementation; this is also where I acquired the array containing the interpolation coefficients).  
My implementation seems to be very fast (when compared to approaches involving, e.g., `gsl`'s interpolation), especially because it handles non-regular grids without direct performance losses. Interpolation speed is independent of the grid spacings because these only play a role in the initial caching of the interpolation coefficients, after which any interpolation runs more or less with the same speed (~40M interpolations per second on my computer, single thread).

## Usage
### Memory, important!
Note that the default fast version of the interpolator (more info below) needs significant amounts of memory (easily >1GB for medium grid sizes) for the cached values. You can estimate the memory needed with
```
mem = nx * ny * nz * 64 * 8 * 1.0e-9 GB,
```
where `nx`, `ny` and `nz` are the sizes of the grid in `x`, `y` and `z` directions.  
A slower version (~15M interpolations s⁻¹, still faster than most other implementations), without coefficient caching and thus much lower memory needs, is available if compiled with `-DINTERP_NO_CACHE`. _Note: This is not yet implemented, so expect large memory needs._


## Setup
The interpolator needs to be set up before use. Calling any of the interpolation functions (see below) before setup will probably just result in a seg fault because it does not check its status for performance reasons.  

There are multiple ways to prepare the interpolator for use:
1. Have the program generate the data (see section **Data generation**) to interpolate on:  
This only really viable if the initial generation of data is slow (for example numerical integration) and if you can export the data (see section **Data export**).
2. Import data (see section **Data import**) that was previously exported:  
This is the recommended way to use the interpolator.
3. Import data from a file that you wrote to a file yourself (see section **Data file format**)

### Data generation
Data generation is accomplished by calling the method `Interpolation3D::generate_data(std::function<double (double, double, double)> f, const DataGenerationConfig* config,  bool monitor_progress = false)`. The `struct DataGenerationConfig` allows you to configure the grid with its fields.  
_In the following explanations, replace the hash symbol `#` with `x`, `y` or `z` to configure the corresponding direction._
* `n#`:  
Is the number of grid points in direction `#`. Ensure that `n# > 4`, as at least 4 points are needed for tricubic interpolation.
* `#_min` and `#_max`:  
Define respectively the minimum and maximum values that will be plugged into the specific function.
* `#_exp_grid_spacing_parameter`:  
Describes the density of grid points in direction `#`. Values `>0` result in greater density at the `#_min`, and smaller density at `#_max`, values `<0` result in the opposite. A value of `0.0` commands uniform density over the whole range.  
_Note: In order for values other than `0.0` to take effect, you also need to set `#_grid_spacing = Exponential`._

### Data export
WIP

