//=============================================================================
//
// nldiffusion_functions.cpp
// Author: Pablo F. Alcantarilla
// Institution: University d'Auvergne
// Address: Clermont Ferrand, France
// Date: 27/12/2011
// Email: pablofdezalc@gmail.com
//
// KAZE Features Copyright 2012, Pablo F. Alcantarilla
// All Rights Reserved
// See LICENSE for the license information
//=============================================================================

/**
 * @file nldiffusion_functions.cpp
 * @brief Functions for non-linear diffusion applications:
 * 2D Gaussian Derivatives
 * Perona and Malik conductivity equations
 * Perona and Malik evolution
 * @date Dec 27, 2011
 * @author Pablo F. Alcantarilla
 */

#include "precomp.hpp"
#include "nldiffusion_functions.h"


/* ************************************************************************* */

namespace libakaze
{
using namespace libakaze;
#define EXPTAB_SCALE 6
#define EXPTAB_MASK  ((1 << EXPTAB_SCALE) - 1)

#define EXPPOLY_32F_A0 .9670371139572337719125840413672004409288e-2

static const double expTab[] = {
    1.0 * EXPPOLY_32F_A0,
    1.0108892860517004600204097905619 * EXPPOLY_32F_A0,
    1.0218971486541166782344801347833 * EXPPOLY_32F_A0,
    1.0330248790212284225001082839705 * EXPPOLY_32F_A0,
    1.0442737824274138403219664787399 * EXPPOLY_32F_A0,
    1.0556451783605571588083413251529 * EXPPOLY_32F_A0,
    1.0671404006768236181695211209928 * EXPPOLY_32F_A0,
    1.0787607977571197937406800374385 * EXPPOLY_32F_A0,
    1.0905077326652576592070106557607 * EXPPOLY_32F_A0,
    1.1023825833078409435564142094256 * EXPPOLY_32F_A0,
    1.1143867425958925363088129569196 * EXPPOLY_32F_A0,
    1.126521618608241899794798643787 * EXPPOLY_32F_A0,
    1.1387886347566916537038302838415 * EXPPOLY_32F_A0,
    1.151189229952982705817759635202 * EXPPOLY_32F_A0,
    1.1637248587775775138135735990922 * EXPPOLY_32F_A0,
    1.1763969916502812762846457284838 * EXPPOLY_32F_A0,
    1.1892071150027210667174999705605 * EXPPOLY_32F_A0,
    1.2021567314527031420963969574978 * EXPPOLY_32F_A0,
    1.2152473599804688781165202513388 * EXPPOLY_32F_A0,
    1.2284805361068700056940089577928 * EXPPOLY_32F_A0,
    1.2418578120734840485936774687266 * EXPPOLY_32F_A0,
    1.2553807570246910895793906574423 * EXPPOLY_32F_A0,
    1.2690509571917332225544190810323 * EXPPOLY_32F_A0,
    1.2828700160787782807266697810215 * EXPPOLY_32F_A0,
    1.2968395546510096659337541177925 * EXPPOLY_32F_A0,
    1.3109612115247643419229917863308 * EXPPOLY_32F_A0,
    1.3252366431597412946295370954987 * EXPPOLY_32F_A0,
    1.3396675240533030053600306697244 * EXPPOLY_32F_A0,
    1.3542555469368927282980147401407 * EXPPOLY_32F_A0,
    1.3690024229745906119296011329822 * EXPPOLY_32F_A0,
    1.3839098819638319548726595272652 * EXPPOLY_32F_A0,
    1.3989796725383111402095281367152 * EXPPOLY_32F_A0,
    1.4142135623730950488016887242097 * EXPPOLY_32F_A0,
    1.4296133383919700112350657782751 * EXPPOLY_32F_A0,
    1.4451808069770466200370062414717 * EXPPOLY_32F_A0,
    1.4609177941806469886513028903106 * EXPPOLY_32F_A0,
    1.476826145939499311386907480374 * EXPPOLY_32F_A0,
    1.4929077282912648492006435314867 * EXPPOLY_32F_A0,
    1.5091644275934227397660195510332 * EXPPOLY_32F_A0,
    1.5255981507445383068512536895169 * EXPPOLY_32F_A0,
    1.5422108254079408236122918620907 * EXPPOLY_32F_A0,
    1.5590044002378369670337280894749 * EXPPOLY_32F_A0,
    1.5759808451078864864552701601819 * EXPPOLY_32F_A0,
    1.5931421513422668979372486431191 * EXPPOLY_32F_A0,
    1.6104903319492543081795206673574 * EXPPOLY_32F_A0,
    1.628027421857347766848218522014 * EXPPOLY_32F_A0,
    1.6457554781539648445187567247258 * EXPPOLY_32F_A0,
    1.6636765803267364350463364569764 * EXPPOLY_32F_A0,
    1.6817928305074290860622509524664 * EXPPOLY_32F_A0,
    1.7001063537185234695013625734975 * EXPPOLY_32F_A0,
    1.7186192981224779156293443764563 * EXPPOLY_32F_A0,
    1.7373338352737062489942020818722 * EXPPOLY_32F_A0,
    1.7562521603732994831121606193753 * EXPPOLY_32F_A0,
    1.7753764925265212525505592001993 * EXPPOLY_32F_A0,
    1.7947090750031071864277032421278 * EXPPOLY_32F_A0,
    1.8142521755003987562498346003623 * EXPPOLY_32F_A0,
    1.8340080864093424634870831895883 * EXPPOLY_32F_A0,
    1.8539791250833855683924530703377 * EXPPOLY_32F_A0,
    1.8741676341102999013299989499544 * EXPPOLY_32F_A0,
    1.8945759815869656413402186534269 * EXPPOLY_32F_A0,
    1.9152065613971472938726112702958 * EXPPOLY_32F_A0,
    1.9360617934922944505980559045667 * EXPPOLY_32F_A0,
    1.9571441241754002690183222516269 * EXPPOLY_32F_A0,
    1.9784560263879509682582499181312 * EXPPOLY_32F_A0,
};
#define EXPPOLY_32F_A0 .9670371139572337719125840413672004409288e-2
static const double exp_prescale = 1.4426950408889634073599246810019 * (1 << EXPTAB_SCALE);
static const double exp_postscale = 1./(1 << EXPTAB_SCALE);
static const double exp_max_val = 3000.*(1 << EXPTAB_SCALE); // log10(DBL_MAX) < 3000
static void Exp_32f( const float *_x, float *y, int n )
{
    static const float
        A4 = (float)(1.000000000000002438532970795181890933776 / EXPPOLY_32F_A0),
        A3 = (float)(.6931471805521448196800669615864773144641 / EXPPOLY_32F_A0),
        A2 = (float)(.2402265109513301490103372422686535526573 / EXPPOLY_32F_A0),
        A1 = (float)(.5550339366753125211915322047004666939128e-1 / EXPPOLY_32F_A0);

#undef EXPPOLY
#define EXPPOLY(x)  \
    (((((x) + A1)*(x) + A2)*(x) + A3)*(x) + A4)

    int i = 0;
    const Cv32suf* x = (const Cv32suf*)_x;
    Cv32suf buf[4];
    for( ; i <= n - 4; i += 4 )
    {
        double x0 = x[i].f * exp_prescale;
        double x1 = x[i + 1].f * exp_prescale;
        double x2 = x[i + 2].f * exp_prescale;
        double x3 = x[i + 3].f * exp_prescale;
        int val0, val1, val2, val3, t;

        if( ((x[i].i >> 23) & 255) > 127 + 10 )
            x0 = x[i].i < 0 ? -exp_max_val : exp_max_val;

        if( ((x[i+1].i >> 23) & 255) > 127 + 10 )
            x1 = x[i+1].i < 0 ? -exp_max_val : exp_max_val;

        if( ((x[i+2].i >> 23) & 255) > 127 + 10 )
            x2 = x[i+2].i < 0 ? -exp_max_val : exp_max_val;

        if( ((x[i+3].i >> 23) & 255) > 127 + 10 )
            x3 = x[i+3].i < 0 ? -exp_max_val : exp_max_val;

        val0 = cvRound(x0);
        val1 = cvRound(x1);
        val2 = cvRound(x2);
        val3 = cvRound(x3);

        x0 = (x0 - val0)*exp_postscale;
        x1 = (x1 - val1)*exp_postscale;
        x2 = (x2 - val2)*exp_postscale;
        x3 = (x3 - val3)*exp_postscale;

        t = (val0 >> EXPTAB_SCALE) + 127;
        t = !(t & ~255) ? t : t < 0 ? 0 : 255;
        buf[0].i = t << 23;

        t = (val1 >> EXPTAB_SCALE) + 127;
        t = !(t & ~255) ? t : t < 0 ? 0 : 255;
        buf[1].i = t << 23;

        t = (val2 >> EXPTAB_SCALE) + 127;
        t = !(t & ~255) ? t : t < 0 ? 0 : 255;
        buf[2].i = t << 23;

        t = (val3 >> EXPTAB_SCALE) + 127;
        t = !(t & ~255) ? t : t < 0 ? 0 : 255;
        buf[3].i = t << 23;

        x0 = buf[0].f * expTab[val0 & EXPTAB_MASK] * EXPPOLY( x0 );
        x1 = buf[1].f * expTab[val1 & EXPTAB_MASK] * EXPPOLY( x1 );

        y[i] = (float)x0;
        y[i + 1] = (float)x1;

        x2 = buf[2].f * expTab[val2 & EXPTAB_MASK] * EXPPOLY( x2 );
        x3 = buf[3].f * expTab[val3 & EXPTAB_MASK] * EXPPOLY( x3 );

        y[i + 2] = (float)x2;
        y[i + 3] = (float)x3;
    }

    for( ; i < n; i++ )
    {
        double x0 = x[i].f * exp_prescale;
        int val0, t;

        if( ((x[i].i >> 23) & 255) > 127 + 10 )
            x0 = x[i].i < 0 ? -exp_max_val : exp_max_val;

        val0 = cvRound(x0);
        t = (val0 >> EXPTAB_SCALE) + 127;
        t = !(t & ~255) ? t : t < 0 ? 0 : 255;

        buf[0].i = t << 23;
        x0 = (x0 - val0)*exp_postscale;

        y[i] = (float)(buf[0].f * expTab[val0 & EXPTAB_MASK] * EXPPOLY(x0));
    }
}


static void twexp( InputArray _src, OutputArray _dst )
{
    Mat src = _src.getMat();
    int type = src.type(), depth = src.depth(), cn = src.channels();

    _dst.create( src.dims, src.size, type );
    Mat dst = _dst.getMat();

    CV_Assert( depth == CV_32F );

    const Mat* arrays[] = {&src, &dst, 0};
    uchar* ptrs[2];
    NAryMatIterator it(arrays, ptrs);
    int len = (int)(it.size*cn);

    for( size_t i = 0; i < it.nplanes; i++, ++it )
    {
        Exp_32f( (const float*)ptrs[0], (float*)ptrs[1], len );
    }
}
/* ************************************************************************* */
/**
 * @brief This function smoothes an image with a Gaussian kernel
 * @param src Input image
 * @param dst Output image
 * @param ksize_x Kernel size in X-direction (horizontal)
 * @param ksize_y Kernel size in Y-direction (vertical)
 * @param sigma Kernel standard deviation
 */
void gaussian_2D_convolution(const libakaze::Mat& src, libakaze::Mat& dst, int ksize_x, int ksize_y, float sigma) {

    int ksize_x_ = 0, ksize_y_ = 0;

    // Compute an appropriate kernel size according to the specified sigma
    if (sigma > ksize_x || sigma > ksize_y || ksize_x == 0 || ksize_y == 0) {
        ksize_x_ = cvCeil(2.0f*(1.0f + (sigma - 0.8f) / (0.3f)));
        ksize_y_ = ksize_x_;
    }

    // The kernel size must be and odd number
    if ((ksize_x_ % 2) == 0) {
        ksize_x_ += 1;
    }

    if ((ksize_y_ % 2) == 0) {
        ksize_y_ += 1;
    }

    // Perform the Gaussian Smoothing with border replication
    GaussianBlur(src, dst, Size(ksize_x_, ksize_y_), sigma, sigma, BORDER_REPLICATE);
}

/* ************************************************************************* */
/**
 * @brief This function computes image derivatives with Scharr kernel
 * @param src Input image
 * @param dst Output image
 * @param xorder Derivative order in X-direction (horizontal)
 * @param yorder Derivative order in Y-direction (vertical)
 * @note Scharr operator approximates better rotation invariance than
 * other stencils such as Sobel. See Weickert and Scharr,
 * A Scheme for Coherence-Enhancing Diffusion Filtering with Optimized Rotation Invariance,
 * Journal of Visual Communication and Image Representation 2002
 */
void image_derivatives_scharr(const libakaze::Mat& src, libakaze::Mat& dst, int xorder, int yorder) {
    Scharr(src, dst, CV_32F, xorder, yorder, 1.0, 0, BORDER_DEFAULT);
}

/* ************************************************************************* */
/**
 * @brief This function computes the Perona and Malik conductivity coefficient g1
 * g1 = exp(-|dL|^2/k^2)
 * @param Lx First order image derivative in X-direction (horizontal)
 * @param Ly First order image derivative in Y-direction (vertical)
 * @param dst Output image
 * @param k Contrast factor parameter
 */
void pm_g1(InputArray _Lx, InputArray _Ly, OutputArray _dst, float k) {
  _dst.create(_Lx.size(), _Lx.type());
  Mat Lx = _Lx.getMat();
  Mat Ly = _Ly.getMat();
  Mat dst = _dst.getMat();

  Size sz = Lx.size();
  float inv_k = 1.0f / (k*k);
  for (int y = 0; y < sz.height; y++) {

    const float* Lx_row = Lx.ptr<float>(y);
    const float* Ly_row = Ly.ptr<float>(y);
    float* dst_row = dst.ptr<float>(y);

    for (int x = 0; x < sz.width; x++) {
      dst_row[x] = (-inv_k*(Lx_row[x]*Lx_row[x] + Ly_row[x]*Ly_row[x]));
    }
  }

  twexp(dst, dst);
}

/* ************************************************************************* */
/**
 * @brief This function computes the Perona and Malik conductivity coefficient g2
 * g2 = 1 / (1 + dL^2 / k^2)
 * @param Lx First order image derivative in X-direction (horizontal)
 * @param Ly First order image derivative in Y-direction (vertical)
 * @param dst Output image
 * @param k Contrast factor parameter
 */
void pm_g2(InputArray _Lx, InputArray _Ly, OutputArray _dst, float k) {
//    CV_INSTRUMENT_REGION()

    _dst.create(_Lx.size(), _Lx.type());
    Mat Lx = _Lx.getMat();
    Mat Ly = _Ly.getMat();
    Mat dst = _dst.getMat();

    Size sz = Lx.size();
    dst.create(sz, Lx.type());
    float k2inv = 1.0f / (k * k);

    for(int y = 0; y < sz.height; y++) {
        const float *Lx_row = Lx.ptr<float>(y);
        const float *Ly_row = Ly.ptr<float>(y);
        float* dst_row = dst.ptr<float>(y);
        for(int x = 0; x < sz.width; x++) {
            dst_row[x] = 1.0f / (1.0f + ((Lx_row[x] * Lx_row[x] + Ly_row[x] * Ly_row[x]) * k2inv));
        }
    }
}
/* ************************************************************************* */
/**
 * @brief This function computes Weickert conductivity coefficient gw
 * @param Lx First order image derivative in X-direction (horizontal)
 * @param Ly First order image derivative in Y-direction (vertical)
 * @param dst Output image
 * @param k Contrast factor parameter
 * @note For more information check the following paper: J. Weickert
 * Applications of nonlinear diffusion in image processing and computer vision,
 * Proceedings of Algorithmy 2000
 */
void weickert_diffusivity(InputArray _Lx, InputArray _Ly, OutputArray _dst, float k) {
  _dst.create(_Lx.size(), _Lx.type());
  Mat Lx = _Lx.getMat();
  Mat Ly = _Ly.getMat();
  Mat dst = _dst.getMat();

  Size sz = Lx.size();
  float inv_k = 1.0f / (k*k);
  for (int y = 0; y < sz.height; y++) {

    const float* Lx_row = Lx.ptr<float>(y);
    const float* Ly_row = Ly.ptr<float>(y);
    float* dst_row = dst.ptr<float>(y);

    for (int x = 0; x < sz.width; x++) {
      float dL = inv_k*(Lx_row[x]*Lx_row[x] + Ly_row[x]*Ly_row[x]);
      dst_row[x] = -3.315f/(dL*dL*dL*dL);
    }
  }

  twexp(dst, dst);
  dst = 1.0 - dst;
}


/* ************************************************************************* */
/**
* @brief This function computes Charbonnier conductivity coefficient gc
* gc = 1 / sqrt(1 + dL^2 / k^2)
* @param Lx First order image derivative in X-direction (horizontal)
* @param Ly First order image derivative in Y-direction (vertical)
* @param dst Output image
* @param k Contrast factor parameter
* @note For more information check the following paper: J. Weickert
* Applications of nonlinear diffusion in image processing and computer vision,
* Proceedings of Algorithmy 2000
*/
void charbonnier_diffusivity(InputArray _Lx, InputArray _Ly, OutputArray _dst, float k) {
  _dst.create(_Lx.size(), _Lx.type());
  Mat Lx = _Lx.getMat();
  Mat Ly = _Ly.getMat();
  Mat dst = _dst.getMat();

  Size sz = Lx.size();
  float inv_k = 1.0f / (k*k);
  for (int y = 0; y < sz.height; y++) {

    const float* Lx_row = Lx.ptr<float>(y);
    const float* Ly_row = Ly.ptr<float>(y);
    float* dst_row = dst.ptr<float>(y);

    for (int x = 0; x < sz.width; x++) {
      float den = sqrt(1.0f+inv_k*(Lx_row[x]*Lx_row[x] + Ly_row[x]*Ly_row[x]));
      dst_row[x] = 1.0f / den;
    }
  }
}


/* ************************************************************************* */
/**
 * @brief This function computes a good empirical value for the k contrast factor
 * given an input image, the percentile (0-1), the gradient scale and the number of
 * bins in the histogram
 * @param img Input image
 * @param perc Percentile of the image gradient histogram (0-1)
 * @param gscale Scale for computing the image gradient histogram
 * @param nbins Number of histogram bins
 * @param ksize_x Kernel size in X-direction (horizontal) for the Gaussian smoothing kernel
 * @param ksize_y Kernel size in Y-direction (vertical) for the Gaussian smoothing kernel
 * @return k contrast factor
 */
float compute_k_percentile(const libakaze::Mat& img, float perc, float gscale, int nbins, int ksize_x, int ksize_y) {
//    CV_INSTRUMENT_REGION()

    int nbin = 0, nelements = 0, nthreshold = 0, k = 0;
    float kperc = 0.0, modg = 0.0;
    float npoints = 0.0;
    float hmax = 0.0;

    // Create the array for the histogram
    std::vector<int> hist(nbins, 0);

    // Create the matrices
    Mat gaussian = Mat::zeros(img.rows, img.cols, CV_32F);
    Mat Lx = Mat::zeros(img.rows, img.cols, CV_32F);
    Mat Ly = Mat::zeros(img.rows, img.cols, CV_32F);

    // Perform the Gaussian convolution
    gaussian_2D_convolution(img, gaussian, ksize_x, ksize_y, gscale);

    // Compute the Gaussian derivatives Lx and Ly
    Scharr(gaussian, Lx, CV_32F, 1, 0, 1, 0, libakaze::BORDER_DEFAULT);
    Scharr(gaussian, Ly, CV_32F, 0, 1, 1, 0, libakaze::BORDER_DEFAULT);

    // Skip the borders for computing the histogram
    for (int i = 1; i < gaussian.rows - 1; i++) {
        const float *lx = Lx.ptr<float>(i);
        const float *ly = Ly.ptr<float>(i);
        for (int j = 1; j < gaussian.cols - 1; j++) {
            modg = lx[j]*lx[j] + ly[j]*ly[j];

            // Get the maximum
            if (modg > hmax) {
                hmax = modg;
            }
        }
    }
    hmax = sqrt(hmax);
    // Skip the borders for computing the histogram
    for (int i = 1; i < gaussian.rows - 1; i++) {
        const float *lx = Lx.ptr<float>(i);
        const float *ly = Ly.ptr<float>(i);
        for (int j = 1; j < gaussian.cols - 1; j++) {
            modg = lx[j]*lx[j] + ly[j]*ly[j];

            // Find the correspondent bin
            if (modg != 0.0) {
                nbin = (int)floor(nbins*(sqrt(modg) / hmax));

                if (nbin == nbins) {
                    nbin--;
                }

                hist[nbin]++;
                npoints++;
            }
        }
    }

    // Now find the perc of the histogram percentile
    nthreshold = (int)(npoints*perc);

    for (k = 0; nelements < nthreshold && k < nbins; k++) {
        nelements = nelements + hist[k];
    }

    if (nelements < nthreshold)  {
        kperc = 0.03f;
    }
    else {
        kperc = hmax*((float)(k) / (float)nbins);
    }

    return kperc;
}

/* ************************************************************************* */
/**
 * @brief This function computes Scharr image derivatives
 * @param src Input image
 * @param dst Output image
 * @param xorder Derivative order in X-direction (horizontal)
 * @param yorder Derivative order in Y-direction (vertical)
 * @param scale Scale factor for the derivative size
 */
void compute_scharr_derivatives(const libakaze::Mat& src, libakaze::Mat& dst, int xorder, int yorder, int scale) {
    Mat kx, ky;
    compute_derivative_kernels(kx, ky, xorder, yorder, scale);
    sepFilter2D(src, dst, CV_32F, kx, ky);
}

/* ************************************************************************* */
/**
 * @brief Compute derivative kernels for sizes different than 3
 * @param _kx Horizontal kernel ues
 * @param _ky Vertical kernel values
 * @param dx Derivative order in X-direction (horizontal)
 * @param dy Derivative order in Y-direction (vertical)
 * @param scale_ Scale factor or derivative size
 */
void compute_derivative_kernels(libakaze::OutputArray _kx, libakaze::OutputArray _ky, int dx, int dy, int scale) {
//    CV_INSTRUMENT_REGION()

    int ksize = 3 + 2 * (scale - 1);

    // The standard Scharr kernel
    if (scale == 1) {
        getDerivKernels(_kx, _ky, dx, dy, 0, true, CV_32F);
        return;
    }

    _kx.create(ksize, 1, CV_32F, -1, true);
    _ky.create(ksize, 1, CV_32F, -1, true);
    Mat kx = _kx.getMat();
    Mat ky = _ky.getMat();
    std::vector<float> kerI;

    float w = 10.0f / 3.0f;
    float norm = 1.0f / (2.0f*scale*(w + 2.0f));

    for (int k = 0; k < 2; k++) {
        Mat* kernel = k == 0 ? &kx : &ky;
        int order = k == 0 ? dx : dy;
        kerI.assign(ksize, 0.0f);

        if (order == 0) {
            kerI[0] = norm, kerI[ksize / 2] = w*norm, kerI[ksize - 1] = norm;
        }
        else if (order == 1) {
            kerI[0] = -1, kerI[ksize / 2] = 0, kerI[ksize - 1] = 1;
        }

        Mat temp(kernel->rows, kernel->cols, CV_32F, &kerI[0]);
        temp.copyTo(*kernel);
    }
}

class Nld_Step_Scalar_Invoker : public libakaze::ParallelLoopBody
{
public:
    Nld_Step_Scalar_Invoker(libakaze::Mat& Ld, const libakaze::Mat& c, libakaze::Mat& Lstep, float _stepsize)
        : _Ld(&Ld)
        , _c(&c)
        , _Lstep(&Lstep)
        , stepsize(_stepsize)
    {
    }

    virtual ~Nld_Step_Scalar_Invoker()
    {

    }

    void operator()(const libakaze::Range& range) const
    {
        libakaze::Mat& Ld = *_Ld;
        const libakaze::Mat& c = *_c;
        libakaze::Mat& Lstep = *_Lstep;

        for (int i = range.start; i < range.end; i++)
        {
            const float *c_prev  = c.ptr<float>(i - 1);
            const float *c_curr  = c.ptr<float>(i);
            const float *c_next  = c.ptr<float>(i + 1);
            const float *ld_prev = Ld.ptr<float>(i - 1);
            const float *ld_curr = Ld.ptr<float>(i);
            const float *ld_next = Ld.ptr<float>(i + 1);

            float *dst  = Lstep.ptr<float>(i);

            for (int j = 1; j < Lstep.cols - 1; j++)
            {
                float xpos = (c_curr[j]   + c_curr[j+1])*(ld_curr[j+1] - ld_curr[j]);
                float xneg = (c_curr[j-1] + c_curr[j])  *(ld_curr[j]   - ld_curr[j-1]);
                float ypos = (c_curr[j]   + c_next[j])  *(ld_next[j]   - ld_curr[j]);
                float yneg = (c_prev[j]   + c_curr[j])  *(ld_curr[j]   - ld_prev[j]);
                dst[j] = 0.5f*stepsize*(xpos - xneg + ypos - yneg);
            }
        }
    }
private:
    libakaze::Mat * _Ld;
    const libakaze::Mat * _c;
    libakaze::Mat * _Lstep;
    float stepsize;
};

/* ************************************************************************* */
/**
* @brief This function performs a scalar non-linear diffusion step
* @param Ld2 Output image in the evolution
* @param c Conductivity image
* @param Lstep Previous image in the evolution
* @param stepsize The step size in time units
* @note Forward Euler Scheme 3x3 stencil
* The function c is a scalar value that depends on the gradient norm
* dL_by_ds = d(c dL_by_dx)_by_dx + d(c dL_by_dy)_by_dy
*/
void nld_step_scalar(libakaze::Mat& Ld, const libakaze::Mat& c, libakaze::Mat& Lstep, float stepsize) {
//    CV_INSTRUMENT_REGION()

    libakaze::parallel_for_(libakaze::Range(1, Lstep.rows - 1), Nld_Step_Scalar_Invoker(Ld, c, Lstep, stepsize), (double)Ld.total()/(1 << 16));

    float xneg, xpos, yneg, ypos;
    float* dst = Lstep.ptr<float>(0);
    const float* cprv = NULL;
    const float* ccur  = c.ptr<float>(0);
    const float* cnxt  = c.ptr<float>(1);
    const float* ldprv = NULL;
    const float* ldcur = Ld.ptr<float>(0);
    const float* ldnxt = Ld.ptr<float>(1);
    for (int j = 1; j < Lstep.cols - 1; j++) {
        xpos = (ccur[j]   + ccur[j+1]) * (ldcur[j+1] - ldcur[j]);
        xneg = (ccur[j-1] + ccur[j])   * (ldcur[j]   - ldcur[j-1]);
        ypos = (ccur[j]   + cnxt[j])   * (ldnxt[j]   - ldcur[j]);
        dst[j] = 0.5f*stepsize*(xpos - xneg + ypos);
    }

    dst = Lstep.ptr<float>(Lstep.rows - 1);
    ccur = c.ptr<float>(Lstep.rows - 1);
    cprv = c.ptr<float>(Lstep.rows - 2);
    ldcur = Ld.ptr<float>(Lstep.rows - 1);
    ldprv = Ld.ptr<float>(Lstep.rows - 2);

    for (int j = 1; j < Lstep.cols - 1; j++) {
        xpos = (ccur[j] + ccur[j+1]) * (ldcur[j+1] - ldcur[j]);
        xneg = (ccur[j-1] + ccur[j]) * (ldcur[j] - ldcur[j-1]);
        yneg = (cprv[j] + ccur[j])   * (ldcur[j] - ldprv[j]);
        dst[j] = 0.5f*stepsize*(xpos - xneg - yneg);
    }

    ccur = c.ptr<float>(1);
    ldcur = Ld.ptr<float>(1);
    cprv = c.ptr<float>(0);
    ldprv = Ld.ptr<float>(0);

    int r0 = Lstep.cols - 1;
    int r1 = Lstep.cols - 2;

    for (int i = 1; i < Lstep.rows - 1; i++) {
        cnxt = c.ptr<float>(i + 1);
        ldnxt = Ld.ptr<float>(i + 1);
        dst = Lstep.ptr<float>(i);

        xpos = (ccur[0] + ccur[1]) * (ldcur[1] - ldcur[0]);
        ypos = (ccur[0] + cnxt[0]) * (ldnxt[0] - ldcur[0]);
        yneg = (cprv[0] + ccur[0]) * (ldcur[0] - ldprv[0]);
        dst[0] = 0.5f*stepsize*(xpos + ypos - yneg);

        xneg = (ccur[r1] + ccur[r0]) * (ldcur[r0] - ldcur[r1]);
        ypos = (ccur[r0] + cnxt[r0]) * (ldnxt[r0] - ldcur[r0]);
        yneg = (cprv[r0] + ccur[r0]) * (ldcur[r0] - ldprv[r0]);
        dst[r0] = 0.5f*stepsize*(-xneg + ypos - yneg);

        cprv = ccur;
        ccur = cnxt;
        ldprv = ldcur;
        ldcur = ldnxt;
    }
    Ld += Lstep;
}



/* ************************************************************************* */
/**
 * @brief This function checks if a given pixel is a maximum in a local neighbourhood
 * @param img Input image where we will perform the maximum search
 * @param dsize Half size of the neighbourhood
 * @param value Response value at (x,y) position
 * @param row Image row coordinate
 * @param col Image column coordinate
 * @param same_img Flag to indicate if the image value at (x,y) is in the input image
 * @return 1->is maximum, 0->otherwise
 */
bool check_maximum_neighbourhood(const libakaze::Mat& img, int dsize, float value, int row, int col, bool same_img) {

    bool response = true;

    for (int i = row - dsize; i <= row + dsize; i++) {
        for (int j = col - dsize; j <= col + dsize; j++) {
            if (i >= 0 && i < img.rows && j >= 0 && j < img.cols) {
                if (same_img == true) {
                    if (i != row || j != col) {
                        if ((*(img.ptr<float>(i)+j)) > value) {
                            response = false;
                            return response;
                        }
                    }
                }
                else {
                    if ((*(img.ptr<float>(i)+j)) > value) {
                        response = false;
                        return response;
                    }
                }
            }
        }
    }

    return response;
}

}
