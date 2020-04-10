/**
 * @file nldiffusion_functions.h
 * @brief Functions for non-linear diffusion applications:
 * 2D Gaussian Derivatives
 * Perona and Malik conductivity equations
 * Perona and Malik evolution
 * @date Dec 27, 2011
 * @author Pablo F. Alcantarilla
 */

#ifndef ___NLDIFFUSION_FUNCTIONS_H__
#define ___NLDIFFUSION_FUNCTIONS_H__

/* ************************************************************************* */
// Declaration of functions

namespace libakaze
{
// Gaussian 2D convolution
void gaussian_2D_convolution(const libakaze::Mat& src, libakaze::Mat& dst, int ksize_x, int ksize_y, float sigma);

// Diffusivity functions
void pm_g1(InputArray Lx, InputArray Ly, OutputArray dst, float k);
void pm_g2(InputArray Lx, InputArray Ly, OutputArray dst, float k);
void weickert_diffusivity(InputArray Lx, InputArray Ly, OutputArray dst, float k);
void charbonnier_diffusivity(InputArray Lx, InputArray Ly, OutputArray dst, float k);

float compute_k_percentile(const libakaze::Mat& img, float perc, float gscale, int nbins, int ksize_x, int ksize_y);

// Image derivatives
void compute_scharr_derivatives(const libakaze::Mat& src, libakaze::Mat& dst, int xorder, int yorder, int scale);
void compute_derivative_kernels(libakaze::OutputArray _kx, libakaze::OutputArray _ky, int dx, int dy, int scale);
void image_derivatives_scharr(const libakaze::Mat& src, libakaze::Mat& dst, int xorder, int yorder);

// Nonlinear diffusion filtering scalar step
void nld_step_scalar(libakaze::Mat& Ld, const libakaze::Mat& c, libakaze::Mat& Lstep, float stepsize);

// For non-maxima suppresion
bool check_maximum_neighbourhood(const libakaze::Mat& img, int dsize, float value, int row, int col, bool same_img);
}

#endif
