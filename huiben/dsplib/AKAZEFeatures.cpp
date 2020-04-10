/**
 * @file AKAZEFeatures.cpp
 * @brief Main class for detecting and describing binary features in an
 * accelerated nonlinear scale space
 * @date Sep 15, 2013
 * @author Pablo F. Alcantarilla, Jesus Nuevo
 */

#include "precomp.hpp"

#include "AKAZEFeatures.h"
#include "nldiffusion_functions.h"
//#include <iostream>

// Namespaces
namespace libakaze
{

#define  divUp(a,b) (((a)+(b)-1)/(b))
using namespace libakaze;

/* ************************************************************************* */
/**
 * @brief AKAZEFeatures constructor with input options
 * @param options AKAZEFeatures configuration options
 * @note This constructor allocates memory for the nonlinear scale space
 */
AKAZEFeatures::AKAZEFeatures(const AKAZEOptions& options) : options_(options) {

  ncycles_ = 0;
  reordering_ = true;

  if (options_.descriptor_size > 0 && options_.descriptor >= AKAZEOptions::DESCRIPTOR_MLDB_UPRIGHT) {
    generateDescriptorSubsample(descriptorSamples_, descriptorBits_, options_.descriptor_size,
                                options_.descriptor_pattern_size, options_.descriptor_channels);
  }

  Allocate_Memory_Evolution();
}

/* ************************************************************************* */
/**
 * @brief This method allocates the memory for the nonlinear diffusion evolution
 */
void AKAZEFeatures::Allocate_Memory_Evolution(void) 
{
  //CV_INSTRUMENT_REGION()

  float rfactor = 0.0f;
  int level_height = 0, level_width = 0;

  // maximum size of the area for the descriptor computation
  float smax = 0.0;
  if (options_.descriptor == AKAZEOptions::DESCRIPTOR_MLDB_UPRIGHT || options_.descriptor == AKAZEOptions::DESCRIPTOR_MLDB) {
    smax = 10.0f*sqrtf(2.0f);
  }
  else if (options_.descriptor == AKAZEOptions::DESCRIPTOR_KAZE_UPRIGHT || options_.descriptor == AKAZEOptions::DESCRIPTOR_KAZE) {
    smax = 12.0f*sqrtf(2.0f);
  }

  // Allocate the dimension of the matrices for the evolution
  for (int i = 0, power = 1; i <= options_.omax - 1; i++, power *= 2) {
    rfactor = 1.0f / power;
    level_height = (int)(options_.img_height*rfactor);
    level_width = (int)(options_.img_width*rfactor);

    // Smallest possible octave and allow one scale if the image is small
    if ((level_width < 80 || level_height < 40) && i != 0) {
      options_.omax = i;
      break;
    }



    for (int j = 0; j < options_.nsublevels; j++) {
      MEvolution step;
      step.size = Size(level_width, level_height);
      step.esigma = options_.soffset*pow(2.f, (float)(j) / (float)(options_.nsublevels) + i);
      step.sigma_size = cvRound(step.esigma * options_.derivative_factor / power);  // In fact sigma_size only depends on j
      step.etime = 0.5f * (step.esigma * step.esigma);
      step.octave = i;
      step.sublevel = j;
      step.octave_ratio = (float)power;
      step.border = cvRound(smax * step.sigma_size) + 1;

      evolution_.push_back(step);
    }
  }
}

/* ************************************************************************* */
/**
 * @brief Computes kernel size for Gaussian smoothing if the image
 * @param sigma Kernel standard deviation
 * @returns kernel size
 */
static inline int getGaussianKernelSize(float sigma) 
{
  // Compute an appropriate kernel size according to the specified sigma
  int ksize = (int)cvCeil(2.0f*(1.0f + (sigma - 0.8f) / (0.3f)));
  ksize |= 1; // kernel should be odd
  return ksize;
}

/* ************************************************************************* */
/**
* @brief This function computes a scalar non-linear diffusion step
* @param Lt Base image in the evolution
* @param Lf Conductivity image
* @param Lstep Output image that gives the difference between the current
* Ld and the next Ld being evolved
* @param row_begin row where to start
* @param row_end last row to fill exclusive. the range is [row_begin, row_end).
* @note Forward Euler Scheme 3x3 stencil
* The function c is a scalar value that depends on the gradient norm
* dL_by_ds = d(c dL_by_dx)_by_dx + d(c dL_by_dy)_by_dy
*/
static inline void
nld_step_scalar_one_lane(const Mat& Lt, const Mat& Lf, Mat& Lstep, float step_size, int row_begin, int row_end)
{
//  CV_INSTRUMENT_REGION()
  /* The labeling scheme for this five star stencil:
   [    a    ]
   [ -1 c +1 ]
   [    b    ]
   */

  Lstep.create(Lt.size(), Lt.type());
  const int cols = Lt.cols - 2;
  int row = row_begin;

  const float *lt_a, *lt_c, *lt_b;
  const float *lf_a, *lf_c, *lf_b;
  float *dst;
  float step_r = 0.f;

  // Process the top row
  if (row == 0) {
    lt_c = Lt.ptr<float>(0) + 1;  /* Skip the left-most column by +1 */
    lf_c = Lf.ptr<float>(0) + 1;
    lt_b = Lt.ptr<float>(1) + 1;
    lf_b = Lf.ptr<float>(1) + 1;

    // fill the corner to prevent uninitialized values
    dst = Lstep.ptr<float>(0);
    dst[0] = 0.0f;
    ++dst;

    for (int j = 0; j < cols; j++) {
      step_r = (lf_c[j] + lf_c[j + 1])*(lt_c[j + 1] - lt_c[j]) +
               (lf_c[j] + lf_c[j - 1])*(lt_c[j - 1] - lt_c[j]) +
               (lf_c[j] + lf_b[j    ])*(lt_b[j    ] - lt_c[j]);
      dst[j] = step_r * step_size;
    }

    // fill the corner to prevent uninitialized values
    dst[cols] = 0.0f;
    ++row;
  }

  // Process the middle rows
  int middle_end = std::min(Lt.rows - 1, row_end);
  for (; row < middle_end; ++row)
  {
    lt_a = Lt.ptr<float>(row - 1);
    lf_a = Lf.ptr<float>(row - 1);
    lt_c = Lt.ptr<float>(row    );
    lf_c = Lf.ptr<float>(row    );
    lt_b = Lt.ptr<float>(row + 1);
    lf_b = Lf.ptr<float>(row + 1);
    dst = Lstep.ptr<float>(row);

    // The left-most column
    step_r = (lf_c[0] + lf_c[1])*(lt_c[1] - lt_c[0]) +
             (lf_c[0] + lf_b[0])*(lt_b[0] - lt_c[0]) +
             (lf_c[0] + lf_a[0])*(lt_a[0] - lt_c[0]);
    dst[0] = step_r * step_size;

    lt_a++; lt_c++; lt_b++;
    lf_a++; lf_c++; lf_b++;
    dst++;

    // The middle columns
    for (int j = 0; j < cols; j++)
    {
      step_r = (lf_c[j] + lf_c[j + 1])*(lt_c[j + 1] - lt_c[j]) +
               (lf_c[j] + lf_c[j - 1])*(lt_c[j - 1] - lt_c[j]) +
               (lf_c[j] + lf_b[j    ])*(lt_b[j    ] - lt_c[j]) +
               (lf_c[j] + lf_a[j    ])*(lt_a[j    ] - lt_c[j]);
      dst[j] = step_r * step_size;
    }

    // The right-most column
    step_r = (lf_c[cols] + lf_c[cols - 1])*(lt_c[cols - 1] - lt_c[cols]) +
             (lf_c[cols] + lf_b[cols    ])*(lt_b[cols    ] - lt_c[cols]) +
             (lf_c[cols] + lf_a[cols    ])*(lt_a[cols    ] - lt_c[cols]);
    dst[cols] = step_r * step_size;
  }

  // Process the bottom row (row == Lt.rows - 1)
  if (row_end == Lt.rows) {
    lt_a = Lt.ptr<float>(row - 1) + 1;  /* Skip the left-most column by +1 */
    lf_a = Lf.ptr<float>(row - 1) + 1;
    lt_c = Lt.ptr<float>(row    ) + 1;
    lf_c = Lf.ptr<float>(row    ) + 1;

    // fill the corner to prevent uninitialized values
    dst = Lstep.ptr<float>(row);
    dst[0] = 0.0f;
    ++dst;

    for (int j = 0; j < cols; j++) {
      step_r = (lf_c[j] + lf_c[j + 1])*(lt_c[j + 1] - lt_c[j]) +
               (lf_c[j] + lf_c[j - 1])*(lt_c[j - 1] - lt_c[j]) +
               (lf_c[j] + lf_a[j    ])*(lt_a[j    ] - lt_c[j]);
      dst[j] = step_r * step_size;
    }

    // fill the corner to prevent uninitialized values
    dst[cols] = 0.0f;
  }
}

class NonLinearScalarDiffusionStep : public ParallelLoopBody
{
public:
  NonLinearScalarDiffusionStep(const Mat& Lt, const Mat& Lf, Mat& Lstep, float step_size)
    : Lt_(&Lt), Lf_(&Lf), Lstep_(&Lstep), step_size_(step_size)
  {}

  void operator()(const Range& range) const
  {
    nld_step_scalar_one_lane(*Lt_, *Lf_, *Lstep_, step_size_, range.start, range.end);
  }

private:
  const Mat* Lt_;
  const Mat* Lf_;
  Mat* Lstep_;
  float step_size_;
};

#ifdef HAVE_OPENCL
static inline bool
ocl_non_linear_diffusion_step(InputArray Lt_, InputArray Lf_, OutputArray Lstep_, float step_size)
{
  if(!Lt_.isContinuous())
    return false;

  UMat Lt = Lt_.getUMat();
  UMat Lf = Lf_.getUMat();
  UMat Lstep = Lstep_.getUMat();

  size_t globalSize[] = {(size_t)Lt.cols, (size_t)Lt.rows};

  ocl::Kernel ker("AKAZE_nld_step_scalar", ocl::features2d::akaze_oclsrc);
  if( ker.empty() )
    return false;

  return ker.args(
    ocl::KernelArg::ReadOnly(Lt),
    ocl::KernelArg::PtrReadOnly(Lf),
    ocl::KernelArg::PtrWriteOnly(Lstep),
    step_size).run(2, globalSize, 0, true);
}
#endif // HAVE_OPENCL

static inline void
non_linear_diffusion_step(InputArray Lt_, InputArray Lf_, OutputArray Lstep_, float step_size)
{
//  CV_INSTRUMENT_REGION()

  Lstep_.create(Lt_.size(), Lt_.type());

  //CV_OCL_RUN(Lt_.isUMat() && Lf_.isUMat() && Lstep_.isUMat(),
//ocl_non_linear_diffusion_step(Lt_, Lf_, Lstep_, step_size));

  Mat Lt = Lt_.getMat();
  Mat Lf = Lf_.getMat();
  Mat Lstep = Lstep_.getMat();
  parallel_for_(Range(0, Lt.rows), NonLinearScalarDiffusionStep(Lt, Lf, Lstep, step_size));
}

/**
 * @brief This function computes a good empirical value for the k contrast factor
 * given two gradient images, the percentile (0-1), the temporal storage to hold
 * gradient norms and the histogram bins
 * @param Lx Horizontal gradient of the input image
 * @param Ly Vertical gradient of the input image
 * @param nbins Number of histogram bins
 * @return k contrast factor
 */
static inline float
compute_kcontrast(InputArray Lx_, InputArray Ly_, float perc, int nbins)
{
  //CV_INSTRUMENT_REGION()

  CV_Assert(nbins > 2);
  CV_Assert(!Lx_.empty());

  Mat Lx = Lx_.getMat();
  Mat Ly = Ly_.getMat();

  // temporary square roots of dot product
  Mat modgs (Lx.rows - 2, Lx.cols - 2, CV_32F);
  const int total = modgs.cols * modgs.rows;
  float *modg = modgs.ptr<float>();
  float hmax = 0.0f;

  for (int i = 1; i < Lx.rows - 1; i++) {
    const float *lx = Lx.ptr<float>(i) + 1;
    const float *ly = Ly.ptr<float>(i) + 1;
    const int cols = Lx.cols - 2;

    for (int j = 0; j < cols; j++) {
      float dist = sqrtf(lx[j] * lx[j] + ly[j] * ly[j]);
      *modg++ = dist;
      hmax = std::max(hmax, dist);
    }
  }
  modg = modgs.ptr<float>();

  if (hmax == 0.0f)
    return 0.03f;  // e.g. a blank image

  // Compute the bin numbers: the value range [0, hmax] -> [0, nbins-1]
  modgs *= (nbins - 1) / hmax;

  // Count up histogram
  std::vector<int> hist(nbins, 0);
  for (int i = 0; i < total; i++)
    hist[(int)modg[i]]++;

  // Now find the perc of the histogram percentile
  const int nthreshold = (int)((total - hist[0]) * perc);  // Exclude hist[0] as background
  int nelements = 0;
  for (int k = 1; k < nbins; k++) {
    if (nelements >= nthreshold)
        return (float)hmax * k / nbins;

    nelements += hist[k];
  }

  return 0.03f;
}


/**
 * @brief Converts input image to grayscale float image
 *
 * @param image any image
 * @param dst grayscale float image
 */
static inline void prepareInputImage(InputArray image, OutputArray dst)
{
  Mat img = image.getMat();
  if (img.channels()== 1)
  {
	   if ( img.depth() == CV_8U )
    img.convertTo(dst, CV_32F, 1.0 / 255.0, 0);
  else if ( img.depth() == CV_16U )
    img.convertTo(dst, CV_32F, 1.0 / 65535.0, 0);
  }
  //  cvtColor(image, img, COLOR_BGR2GRAY);


}

/**
 * @brief This method creates the nonlinear scale space for a given image
 * @param image Input image for which the nonlinear scale space needs to be created
 */
template<typename MatType>
static inline void
create_nonlinear_scale_space(InputArray image, const AKAZEOptions &options,
  const std::vector<std::vector<float > > &tsteps_evolution, std::vector<Evolution<MatType> > &evolution)
{
  //CV_INSTRUMENT_REGION()
  //CV_Assert(evolution.size() > 0);

  // convert input to grayscale float image if needed
  MatType img;
  prepareInputImage(image, img);

  // create first level of the evolution
  int ksize = getGaussianKernelSize(options.soffset);
  GaussianBlur(img, evolution[0].Lsmooth, Size(ksize, ksize), options.soffset, options.soffset, BORDER_REPLICATE);
  evolution[0].Lsmooth.copyTo(evolution[0].Lt);

  if (evolution.size() == 1) 
  {
    // we don't need to compute kcontrast factor
    Compute_Determinant_Hessian_Response(evolution);
    return;
  }

  return;
}

/**
 * @brief Converts between UMatPyramid and Pyramid and vice versa
 * @details Matrices in evolution levels will be copied
 *
 * @param src source pyramid
 * @param dst destination pyramid
 */
template<typename MatTypeSrc, typename MatTypeDst>
static inline void
convertScalePyramid(const std::vector<Evolution<MatTypeSrc> >& src, std::vector<Evolution<MatTypeDst> > &dst)
{
  dst.resize(src.size());
  for (size_t i = 0; i < src.size(); ++i) {
    dst[i] = Evolution<MatTypeDst>(src[i]);
  }
}

/**
 * @brief This method creates the nonlinear scale space for a given image
 * @param image Input image for which the nonlinear scale space needs to be created
 */
void AKAZEFeatures::Create_Nonlinear_Scale_Space(InputArray image)
{
    create_nonlinear_scale_space(image, options_, tsteps_, evolution_);
}

/* ************************************************************************* */



/**
 * @brief Compute determinant from hessians
 * @details Compute Ldet by (Lxx.mul(Lyy) - Lxy.mul(Lxy)) * sigma
 *
 * @param Lxx spatial derivates
 * @param Lxy spatial derivates
 * @param Lyy spatial derivates
 * @param Ldet output determinant
 * @param sigma determinant will be scaled by this sigma
 */
static inline void compute_determinant(InputArray Lxx_, InputArray Lxy_, InputArray Lyy_,
  OutputArray Ldet_, float sigma)
{
//  CV_INSTRUMENT_REGION()

  Ldet_.create(Lxx_.size(), Lxx_.type());

  //CV_OCL_RUN(Lxx_.isUMat() && Ldet_.isUMat(), ocl_compute_determinant(Lxx_, Lxy_, Lyy_, Ldet_, sigma));

  // output determinant
  Mat Lxx = Lxx_.getMat(), Lxy = Lxy_.getMat(), Lyy = Lyy_.getMat(), Ldet = Ldet_.getMat();
  float *lxx = Lxx.ptr<float>();
  float *lxy = Lxy.ptr<float>();
  float *lyy = Lyy.ptr<float>();
  float *ldet = Ldet.ptr<float>();
  const int total = Lxx.cols * Lxx.rows;
  for (int j = 0; j < total; j++) {
    ldet[j] = (lxx[j] * lyy[j] - lxy[j] * lxy[j]) * sigma;
  }

}

template <typename MatType>
class DeterminantHessianResponse : public ParallelLoopBody
{
public:
    explicit DeterminantHessianResponse(std::vector<Evolution<MatType> >& ev)
    : evolution_(&ev)
  {
  }

  void operator()(const Range& range) const
  {
    MatType Lxx, Lxy, Lyy;

    for (int i = range.start; i < range.end; i++)
    {
      Evolution<MatType> &e = (*evolution_)[i];

      // we cannot use cv:Scharr here, because we need to handle also
      // kernel sizes other than 3, by default we are using 9x9, 5x5 and 7x7

      // compute kernels
      Mat DxKx, DxKy, DyKx, DyKy;
      compute_derivative_kernels(DxKx, DxKy, 1, 0, e.sigma_size);
      compute_derivative_kernels(DyKx, DyKy, 0, 1, e.sigma_size);

      // compute the multiscale derivatives
      sepFilter2D(e.Lsmooth, e.Lx, CV_32F, DxKx, DxKy);
      sepFilter2D(e.Lx, Lxx, CV_32F, DxKx, DxKy);
      sepFilter2D(e.Lx, Lxy, CV_32F, DyKx, DyKy);
      sepFilter2D(e.Lsmooth, e.Ly, CV_32F, DyKx, DyKy);
      sepFilter2D(e.Ly, Lyy, CV_32F, DyKx, DyKy);

      // free Lsmooth to same some space in the pyramid, it is not needed anymore
      e.Lsmooth.release();

      // compute determinant scaled by sigma
      float sigma_size_quat = (float)(e.sigma_size * e.sigma_size * e.sigma_size * e.sigma_size);
      compute_determinant(Lxx, Lxy, Lyy, e.Ldet, sigma_size_quat);
    }
  }

private:
  std::vector<Evolution<MatType> >*  evolution_;
};


/**
 * @brief This method computes the feature detector response for the nonlinear scale space
 * @details OCL version
 * @note We use the Hessian determinant as the feature detector response

static inline void
Compute_Determinant_Hessian_Response(UMatPyramid &evolution) {
  CV_INSTRUMENT_REGION()

  DeterminantHessianResponse<UMat> body (evolution);
  body(Range(0, (int)evolution.size()));
}
*/
/**
 * @brief This method computes the feature detector response for the nonlinear scale space
 * @details CPU version
 * @note We use the Hessian determinant as the feature detector response
 */
static inline void
Compute_Determinant_Hessian_Response(Pyramid &evolution) {

  parallel_for_(Range(0, (int)evolution.size()), DeterminantHessianResponse<Mat>(evolution));
}

/* ************************************************************************* */

/**
 * @brief This method selects interesting keypoints through the nonlinear scale space
 * @param kpts Vector of detected keypoints
 */
void AKAZEFeatures::Feature_Detection(std::vector<KeyPoint>& kpts)
{
  //CV_INSTRUMENT_REGION()

  kpts.clear();
  std::vector<Mat> keypoints_by_layers;
  Find_Scale_Space_Extrema(keypoints_by_layers);
  Do_Subpixel_Refinement(keypoints_by_layers, kpts);
  Compute_Keypoints_Orientation(kpts);
}

/**
 * @brief This method searches v for a neighbor point of the point candidate p
 * @param x Coordinates of the keypoint candidate to search a neighbor
 * @param y Coordinates of the keypoint candidate to search a neighbor
 * @param mask Matrix holding keypoints positions
 * @param search_radius neighbour radius for searching keypoints
 * @param idx The index to mask, pointing to keypoint found.
 * @return true if a neighbor point is found; false otherwise
 */
static inline bool
find_neighbor_point(const int x, const int y, const Mat &mask, const int search_radius, int &idx)
{
  // search neighborhood for keypoints
  for (int i = y - search_radius; i < y + search_radius; ++i) {
    const uchar *curr = mask.ptr<uchar>(i);
    for (int j = x - search_radius; j < x + search_radius; ++j) {
      if (curr[j] == 0) {
        continue; // skip non-keypoint
      }
      // fine-compare with L2 metric (L2 is smaller than our search window)
      int dx = j - x;
      int dy = i - y;
      if (dx * dx + dy * dy <= search_radius * search_radius) {
        idx = i * mask.cols + j;
        return true;
      }
    }
  }

  return false;
}

/**
 * @brief Find keypoints in parallel for each pyramid layer
 */
class FindKeypointsSameScale : public ParallelLoopBody
{
public:
    explicit FindKeypointsSameScale(const Pyramid& ev,
      std::vector<Mat>& kpts, float dthreshold)
    : evolution_(&ev), keypoints_by_layers_(&kpts), dthreshold_(dthreshold)
  {}

  void operator()(const Range& range) const
  {
    for (int i = range.start; i < range.end; i++)
    {
      const MEvolution &e = (*evolution_)[i];
      Mat &kpts = (*keypoints_by_layers_)[i];
      // this mask will hold positions of keypoints in this level
      kpts = Mat::zeros(e.Ldet.size(), CV_8UC1);

      // if border is too big we shouldn't search any keypoints
      if (e.border + 1 >= e.Ldet.rows)
        continue;

      const float * prev = e.Ldet.ptr<float>(e.border - 1);
      const float * curr = e.Ldet.ptr<float>(e.border    );
      const float * next = e.Ldet.ptr<float>(e.border + 1);
      const float * ldet = e.Ldet.ptr<float>();
      uchar *mask = kpts.ptr<uchar>();
      const int search_radius = e.sigma_size; // size of keypoint in this level

      for (int y = e.border; y < e.Ldet.rows - e.border; y++) {
        for (int x = e.border; x < e.Ldet.cols - e.border; x++) {
          const float value = curr[x];

          // Filter the points with the detector threshold
          if (value <= dthreshold_)
            continue;
          if (value <= curr[x-1] || value <= curr[x+1])
            continue;
          if (value <= prev[x-1] || value <= prev[x  ] || value <= prev[x+1])
            continue;
          if (value <= next[x-1] || value <= next[x  ] || value <= next[x+1])
            continue;

          int idx = 0;
          // Compare response with the same scale
          if (find_neighbor_point(x, y, kpts, search_radius, idx)) {
            if (value > ldet[idx]) {
              mask[idx] = 0;  // clear old point - we have better candidate now
            } else {
              continue; // there already is a better keypoint
            }
          }

          kpts.at<uchar>(y, x) = 1; // we have a new keypoint
        }

        prev = curr;
        curr = next;
        next += e.Ldet.cols;
      }
    }
  }

private:
  const Pyramid*  evolution_;
  std::vector<Mat>* keypoints_by_layers_;
  float dthreshold_; ///< Detector response threshold to accept point
};

/**
 * @brief This method finds extrema in the nonlinear scale space
 * @param keypoints_by_layers Output vectors of detected keypoints; one vector for each evolution level
 */
void AKAZEFeatures::Find_Scale_Space_Extrema(std::vector<Mat>& keypoints_by_layers)
{
  //CV_INSTRUMENT_REGION()

  keypoints_by_layers.resize(evolution_.size());

  // find points in the same level
  parallel_for_(Range(0, (int)evolution_.size()),FindKeypointsSameScale(evolution_, keypoints_by_layers, options_.dthreshold));
}

/* ************************************************************************* */



/**
 * @brief This method performs subpixel refinement of the detected keypoints
 * @param keypoints_by_layers Input vectors of detected keypoints, sorted by evolution levels
 * @param kpts Output vector of the final refined keypoints
 */
static inline bool tw_solve(float*_src, float*_src2arg, float*_dst)
{
 
	bool result = true;
	#define Sf( y, x ) _src[y*2+x]
	#define bf(y) _src2arg[y]
	#define Df( y, x ) _dst[y]
	#define det2(m)   ((double)m(0,0)*m(1,1) - (double)m(0,1)*m(1,0))
	double d = det2(Sf);
	if( d != 0. )
	{
		double t;
		d = 1./d;
		t = (float)(((double)bf(0)*Sf(1,1) - (double)bf(1)*Sf(0,1))*d);
		Df(1,0) = (float)(((double)bf(1)*Sf(0,0) - (double)bf(0)*Sf(1,0))*d);
		Df(0,0) = (float)t;
	}
	else
		result = false;
   return result;     
}
void AKAZEFeatures::Do_Subpixel_Refinement(
  std::vector<Mat>& keypoints_by_layers, std::vector<KeyPoint>& output_keypoints)
{
  //CV_INSTRUMENT_REGION()

  for (size_t i = 0; i < keypoints_by_layers.size(); i++) {
    const MEvolution &e = evolution_[i];
    const float * const ldet = e.Ldet.ptr<float>();
    const float ratio = e.octave_ratio;
    const int cols = e.Ldet.cols;
    const Mat& keypoints = keypoints_by_layers[i];
    const uchar *const kpts = keypoints.ptr<uchar>();

    size_t j = 0;
    for (int y = 0; y < keypoints.rows; y++) {
      for (int x = 0; x < keypoints.cols; x++, j++) {
        if (kpts[j] == 0) {
          continue; // skip non-keypoints
        }

        // create a new keypoint
        KeyPoint kp;
        kp.pt.x = x * e.octave_ratio;
        kp.pt.y = y * e.octave_ratio;
        kp.size = e.esigma * options_.derivative_factor;
        kp.angle = -1;
        kp.response = ldet[j];
        kp.octave = e.octave;
        kp.class_id = static_cast<int>(i);

        // Compute the gradient
        float Dx = 0.5f * (ldet[ y     *cols + x + 1] - ldet[ y     *cols + x - 1]);
        float Dy = 0.5f * (ldet[(y + 1)*cols + x    ] - ldet[(y - 1)*cols + x    ]);

        // Compute the Hessian
        float Dxx = ldet[ y     *cols + x + 1] + ldet[ y     *cols + x - 1] - 2.0f * ldet[y*cols + x];
        float Dyy = ldet[(y + 1)*cols + x    ] + ldet[(y - 1)*cols + x    ] - 2.0f * ldet[y*cols + x];
        float Dxy = 0.25f * (ldet[(y + 1)*cols + x + 1] + ldet[(y - 1)*cols + x - 1] -
                            ldet[(y - 1)*cols + x + 1] - ldet[(y + 1)*cols + x - 1]);

        // Solve the linear system
      //  Matx22f A( Dxx, Dxy,
       //            Dxy, Dyy );
      //  Vec2f   b( -Dx, -Dy );
      //  Vec2f   dst( 0.0f, 0.0f );
	float _src[2][2]={Dxx,Dxy,Dxy,Dyy};
	float _dst[2],_src2arg[2] = {-Dx,-Dy};
   // solve(A, b, dst, DECOMP_LU);
	tw_solve(_src[0], _src2arg, _dst);
	//printf("%f=,%f=\n",10000000*(_dst[0] - dst[0]),10000000*(_dst[1] - dst[1]));
        float dx = _dst[0];
        float dy = _dst[1];

        if (fabs(dx) > 1.0f || fabs(dy) > 1.0f)
          continue; // Ignore the point that is not stable

        // Refine the coordinates
        kp.pt.x += dx * ratio + .5f*(ratio-1.f);
        kp.pt.y += dy * ratio + .5f*(ratio-1.f);

        kp.angle = 0.0;
        kp.size *= 2.0f; // In OpenCV the size of a keypoint is the diameter

        // Push the refined keypoint to the final storage
        output_keypoints.push_back(kp);
      }
    }
  }
}

/* ************************************************************************* */


class MLDB_Descriptor_Subset_Invoker : public ParallelLoopBody
{
public:
  MLDB_Descriptor_Subset_Invoker(std::vector<KeyPoint>& kpts,
                                 Mat& desc,
                                 Pyramid& evolution,
                                 AKAZEOptions& options,
                                 Mat descriptorSamples,
                                 Mat descriptorBits)
    : keypoints_(&kpts)
    , descriptors_(&desc)
    , evolution_(&evolution)
    , options_(&options)
    , descriptorSamples_(descriptorSamples)
    , descriptorBits_(descriptorBits)
  {
  }

  void operator() (const Range& range) const
  {
    for (int i = range.start; i < range.end; i++)
    {
      Get_MLDB_Descriptor_Subset((*keypoints_)[i], descriptors_->ptr<unsigned char>(i), descriptors_->cols);
    }
  }

  void Get_MLDB_Descriptor_Subset(const KeyPoint& kpt, unsigned char* desc, int desc_size) const;

private:
  std::vector<KeyPoint>* keypoints_;
  Mat*                   descriptors_;
  Pyramid*   evolution_;
  AKAZEOptions*              options_;

  Mat descriptorSamples_;  // List of positions in the grids to sample LDB bits from.
  Mat descriptorBits_;
};

/**
 * @brief This method  computes the set of descriptors through the nonlinear scale space
 * @param kpts Vector of detected keypoints
 * @param desc Matrix to store the descriptors
 */
void AKAZEFeatures::Compute_Descriptors(std::vector<KeyPoint>& kpts, OutputArray descriptors)
{
  //CV_INSTRUMENT_REGION()
  // Allocate memory for the matrix with the descriptors
  int descriptor_size = 32;
  int descriptor_type = CV_8UC1;
  descriptors.create((int)kpts.size(), descriptor_size, descriptor_type);
  Mat desc = descriptors.getMat();
  parallel_for_(Range(0, (int)kpts.size()), MLDB_Descriptor_Subset_Invoker(kpts, desc, evolution_, options_, descriptorSamples_, descriptorBits_));

}

/* ************************************************************************* */
/**
 * @brief This function samples the derivative responses Lx and Ly for the points
 * within the radius of 6*scale from (x0, y0), then multiply 2D Gaussian weight
 * @param Lx Horizontal derivative
 * @param Ly Vertical derivative
 * @param x0 X-coordinate of the center point
 * @param y0 Y-coordinate of the center point
 * @param scale The sampling step
 * @param resX Output array of the weighted horizontal derivative responses
 * @param resY Output array of the weighted vertical derivative responses
 */
static inline
void Sample_Derivative_Response_Radius6(const Mat &Lx, const Mat &Ly,
                                  const int x0, const int y0, const int scale,
                                  float *resX, float *resY)
{
    /* ************************************************************************* */
    /// Lookup table for 2d gaussian (sigma = 2.5) where (0,0) is top left and (6,6) is bottom right
    static const float gauss25[7][7] =
    {
        { 0.02546481f, 0.02350698f, 0.01849125f, 0.01239505f, 0.00708017f, 0.00344629f, 0.00142946f },
        { 0.02350698f, 0.02169968f, 0.01706957f, 0.01144208f, 0.00653582f, 0.00318132f, 0.00131956f },
        { 0.01849125f, 0.01706957f, 0.01342740f, 0.00900066f, 0.00514126f, 0.00250252f, 0.00103800f },
        { 0.01239505f, 0.01144208f, 0.00900066f, 0.00603332f, 0.00344629f, 0.00167749f, 0.00069579f },
        { 0.00708017f, 0.00653582f, 0.00514126f, 0.00344629f, 0.00196855f, 0.00095820f, 0.00039744f },
        { 0.00344629f, 0.00318132f, 0.00250252f, 0.00167749f, 0.00095820f, 0.00046640f, 0.00019346f },
        { 0.00142946f, 0.00131956f, 0.00103800f, 0.00069579f, 0.00039744f, 0.00019346f, 0.00008024f }
    };
    static const struct gtable
    {
      float weight[109];
      int xidx[109];
      int yidx[109];

      explicit gtable(void)
      {
        // Generate the weight and indices by one-time initialization
        int k = 0;
        for (int i = -6; i <= 6; ++i) {
          for (int j = -6; j <= 6; ++j) {
            if (i*i + j*j < 36) {
              CV_Assert(k < 109);
              weight[k] = gauss25[std::abs(i)][std::abs(j)];
              yidx[k] = i;
              xidx[k] = j;
              ++k;
            }
          }
        }
      }
    } g;

  CV_Assert(x0 - 6 * scale >= 0 && x0 + 6 * scale < Lx.cols);
  CV_Assert(y0 - 6 * scale >= 0 && y0 + 6 * scale < Lx.rows);

  for (int i = 0; i < 109; i++)
  {
    int y = y0 + g.yidx[i] * scale;
    int x = x0 + g.xidx[i] * scale;

    float w = g.weight[i];
    resX[i] = w * Lx.at<float>(y, x);
    resY[i] = w * Ly.at<float>(y, x);
  }
}

/**
 * @brief This function sorts a[] by quantized float values
 * @param a[] Input floating point array to sort
 * @param n The length of a[]
 * @param quantum The interval to convert a[i]'s float values to integers
 * @param nkeys a[i] < nkeys * quantum
 * @param idx[] Output array of the indices: a[idx[i]] forms a sorted array
 * @param cum[] Output array of the starting indices of quantized floats
 * @note The values of a[] in [k*quantum, (k + 1)*quantum) is labeled by
 * the integer k, which is calculated by floor(a[i]/quantum).  After sorting,
 * the values from a[idx[cum[k]]] to a[idx[cum[k+1]-1]] are all labeled by k.
 * This sorting is unstable to reduce the memory access.
 */
static const float atan2_p1 = 0.9997878412794807f*(float)(180 / CV_PI);
static const float atan2_p3 = -0.3258083974640975f*(float)(180 / CV_PI);
static const float atan2_p5 = 0.1555786518463281f*(float)(180 / CV_PI);
static const float atan2_p7 = -0.04432655554792128f*(float)(180 / CV_PI);
static inline void FastAtan2_32f(const float *Y, const float *X, float *angle, int len, bool angleInDegrees = true)
{
	int i = 0;
	float scale = angleInDegrees ? 1 : (float)(CV_PI / 180);

	for (; i < len; i++)
	{
		float x = X[i], y = Y[i];
		float ax = std::abs(x), ay = std::abs(y);
		float a, c, c2;
		if (ax >= ay)
		{
			c = ay / (ax + (float)DBL_EPSILON);
			c2 = c*c;
			a = (((atan2_p7*c2 + atan2_p5)*c2 + atan2_p3)*c2 + atan2_p1)*c;
		}
		else
		{
			c = ax / (ay + (float)DBL_EPSILON);
			c2 = c*c;
			a = 90.f - (((atan2_p7*c2 + atan2_p5)*c2 + atan2_p3)*c2 + atan2_p1)*c;
		}
		if (x < 0)
			a = 180.f - a;
		if (y < 0)
			a = 360.f - a;
		angle[i] = (float)(a*scale);
	}
}
static inline
void quantized_counting_sort(const float a[], const int n,
                             const float quantum, const int nkeys,
                             int idx[/*n*/], int cum[/*nkeys + 1*/])
{
  memset(cum, 0, sizeof(cum[0]) * (nkeys + 1));

  // Count up the quantized values
  for (int i = 0; i < n; i++)
  {
    int b = (int)(a[i] / quantum);
    if (b < 0 || b >= nkeys)
      b = 0;
    cum[b]++;
  }

  // Compute the inclusive prefix sum i.e. the end indices; cum[nkeys] is the total
  for (int i = 1; i <= nkeys; i++)
  {
    cum[i] += cum[i - 1];
  }
  CV_Assert(cum[nkeys] == n);

  // Generate the sorted indices; cum[] becomes the exclusive prefix sum i.e. the start indices of keys
  for (int i = 0; i < n; i++)
  {
    int b = (int)(a[i] / quantum);
    if (b < 0 || b >= nkeys)
      b = 0;
    idx[--cum[b]] = i;
  }
}

/**
 * @brief This function computes the main orientation for a given keypoint
 * @param kpt Input keypoint
 * @note The orientation is computed using a similar approach as described in the
 * original SURF method. See Bay et al., Speeded Up Robust Features, ECCV 2006
 */
static inline float twfastAtan2( float y, float x )
{
    float ax = std::abs(x), ay = std::abs(y);
    float a, c, c2;
    if( ax >= ay )
    {
        c = ay/(ax + (float)DBL_EPSILON);
        c2 = c*c;
        a = (((atan2_p7*c2 + atan2_p5)*c2 + atan2_p3)*c2 + atan2_p1)*c;
    }
    else
    {
        c = ax/(ay + (float)DBL_EPSILON);
        c2 = c*c;
        a = 90.f - (((atan2_p7*c2 + atan2_p5)*c2 + atan2_p3)*c2 + atan2_p1)*c;
    }
    if( x < 0 )
        a = 180.f - a;
    if( y < 0 )
        a = 360.f - a;
    return a;
}
static inline
void Compute_Main_Orientation(KeyPoint& kpt, const Pyramid& evolution)
{
  // get the right evolution level for this keypoint
  const MEvolution& e = evolution[kpt.class_id];
  // Get the information from the keypoint
  int scale = cvRound(0.5f * kpt.size / e.octave_ratio);
  int x0 = cvRound(kpt.pt.x / e.octave_ratio);
  int y0 = cvRound(kpt.pt.y / e.octave_ratio);

  // Sample derivatives responses for the points within radius of 6*scale
  const int ang_size = 109;
  float resX[ang_size], resY[ang_size];
  Sample_Derivative_Response_Radius6(e.Lx, e.Ly, x0, y0, scale, resX, resY);

  // Compute the angle of each gradient vector
  float Ang[ang_size];
  //hal::fastAtan2(resY, resX, Ang, ang_size, false);
  FastAtan2_32f(resY, resX, Ang, ang_size, false);

  // Sort by the angles; angles are labeled by slices of 0.15 radian
  const int slices = 42;
  const float ang_step = (float)(2.0 * CV_PI / slices);
  int slice[slices + 1];
  int sorted_idx[ang_size];
  quantized_counting_sort(Ang, ang_size, ang_step, slices, sorted_idx, slice);

  // Find the main angle by sliding a window of 7-slice size(=PI/3) around the keypoint
  const int win = 7;

  float maxX = 0.0f, maxY = 0.0f;
  for (int i = slice[0]; i < slice[win]; i++) {
    const int idx = sorted_idx[i];
    maxX += resX[idx];
    maxY += resY[idx];
  }
  float maxNorm = maxX * maxX + maxY * maxY;

  for (int sn = 1; sn <= slices - win; sn++) {

    if (slice[sn] == slice[sn - 1] && slice[sn + win] == slice[sn + win - 1])
      continue;  // The contents of the window didn't change; don't repeat the computation

    float sumX = 0.0f, sumY = 0.0f;
    for (int i = slice[sn]; i < slice[sn + win]; i++) {
      const int idx = sorted_idx[i];
      sumX += resX[idx];
      sumY += resY[idx];
    }

    float norm = sumX * sumX + sumY * sumY;
    if (norm > maxNorm)
        maxNorm = norm, maxX = sumX, maxY = sumY;  // Found bigger one; update
  }

  for (int sn = slices - win + 1; sn < slices; sn++) {
    int remain = sn + win - slices;

    if (slice[sn] == slice[sn - 1] && slice[remain] == slice[remain - 1])
      continue;

    float sumX = 0.0f, sumY = 0.0f;
    for (int i = slice[sn]; i < slice[slices]; i++) {
      const int idx = sorted_idx[i];
      sumX += resX[idx];
      sumY += resY[idx];
    }
    for (int i = slice[0]; i < slice[remain]; i++) {
      const int idx = sorted_idx[i];
      sumX += resX[idx];
      sumY += resY[idx];
    }

    float norm = sumX * sumX + sumY * sumY;
    if (norm > maxNorm)
        maxNorm = norm, maxX = sumX, maxY = sumY;
  }

  // Store the final result
  kpt.angle = twfastAtan2(maxY, maxX);
}

class ComputeKeypointOrientation : public ParallelLoopBody
{
public:
  ComputeKeypointOrientation(std::vector<KeyPoint>& kpts,
                             const Pyramid& evolution)
    : keypoints_(&kpts)
    , evolution_(&evolution)
  {
  }

  void operator() (const Range& range) const
  {
    for (int i = range.start; i < range.end; i++)
    {
      Compute_Main_Orientation((*keypoints_)[i], *evolution_);
    }
  }
private:
  std::vector<KeyPoint>* keypoints_;
  const Pyramid* evolution_;
};

/**
 * @brief This method computes the main orientation for a given keypoints
 * @param kpts Input keypoints
 */
void AKAZEFeatures::Compute_Keypoints_Orientation(std::vector<KeyPoint>& kpts) const
{
  //CV_INSTRUMENT_REGION()
  parallel_for_(Range(0, (int)kpts.size()), ComputeKeypointOrientation(kpts, evolution_));
}
/* ************************************************************************* */
/**
 * @brief This method computes the M-LDB descriptor of the provided keypoint given the
 * main orientation of the keypoint. The descriptor is computed based on a subset of
 * the bits of the whole descriptor
 * @param kpt Input keypoint
 * @param desc Descriptor vector
 */
void MLDB_Descriptor_Subset_Invoker::Get_MLDB_Descriptor_Subset(const KeyPoint& kpt, unsigned char *desc, int desc_size) const {

  float rx = 0.f, ry = 0.f;
  float sample_x = 0.f, sample_y = 0.f;

  const AKAZEOptions & options = *options_;
  const Pyramid& evolution = *evolution_;

  // Get the information from the keypoint
  float ratio = (float)(1 << kpt.octave);
  int scale = cvRound(0.5f*kpt.size / ratio);
  float angle = kpt.angle * static_cast<float>(CV_PI / 180.f);
  const int level = kpt.class_id;
  const Mat Lx = evolution[level].Lx;
  const Mat Ly = evolution[level].Ly;
  const Mat Lt = evolution[level].Lt;
  float yf = kpt.pt.y / ratio;
  float xf = kpt.pt.x / ratio;
  float co = cos(angle);
  float si = sin(angle);

  // Allocate memory for the matrix of values
  // Buffer for the M-LDB descriptor
  const int max_channels = 3;
  const int channels = options.descriptor_channels;
  CV_Assert(channels <= max_channels);
  float values[(4 + 9 + 16)*max_channels] = { 0 };

  // Sample everything, but only do the comparisons
  const int pattern_size = options.descriptor_pattern_size;
  CV_Assert((pattern_size & 1) == 0);
  const int sample_steps[3] = {
    pattern_size,
    divUp(pattern_size * 2, 3),
    divUp(pattern_size, 2)
  };

  for (int i = 0; i < descriptorSamples_.rows; i++) {
    const int *coords = descriptorSamples_.ptr<int>(i);
    CV_Assert(coords[0] >= 0 && coords[0] < 3);
    const int sample_step = sample_steps[coords[0]];
    float di = 0.f, dx = 0.f, dy = 0.f;

    for (int k = coords[1]; k < coords[1] + sample_step; k++) {
      for (int l = coords[2]; l < coords[2] + sample_step; l++) {

        // Get the coordinates of the sample point
        sample_y = yf + (l*scale*co + k*scale*si);
        sample_x = xf + (-l*scale*si + k*scale*co);

        const int y1 = cvRound(sample_y);
        const int x1 = cvRound(sample_x);

        if (x1 < 0 || y1 < 0 || x1 >= Lt.cols || y1 >= Lt.rows)
            continue; // Boundaries

        di += Lt.at<float>(y1, x1);

        if (options.descriptor_channels > 1) {
          rx = Lx.at<float>(y1, x1);
          ry = Ly.at<float>(y1, x1);

          if (options.descriptor_channels == 2) {
            dx += sqrtf(rx*rx + ry*ry);
          }
          else if (options.descriptor_channels == 3) {
            // Get the x and y derivatives on the rotated axis
            dx += rx*co + ry*si;
            dy += -rx*si + ry*co;
          }
        }
      }
    }

    float* pValues = &values[channels * i];
    pValues[0] = di;

    if (channels == 2) {
      pValues[1] = dx;
    }
    else if (channels == 3) {
      pValues[1] = dx;
      pValues[2] = dy;
    }
  }

  // Do the comparisons
  const int *comps = descriptorBits_.ptr<int>(0);

 // CV_Assert(divUp(descriptorBits_.rows, 8) == desc_size);
  memset(desc, 0, desc_size);

  for (int i = 0; i<descriptorBits_.rows; i++) {
    if (values[comps[2 * i]] > values[comps[2 * i + 1]]) {
      desc[i / 8] |= (1 << (i % 8));
    }
  }
}

/* ************************************************************************* */
/**
 * @brief This method computes the upright (not rotation invariant) M-LDB descriptor
 * of the provided keypoint given the main orientation of the keypoint.
 * The descriptor is computed based on a subset of the bits of the whole descriptor
 * @param kpt Input keypoint
 * @param desc Descriptor vector
 */

/* ************************************************************************* */
/**
 * @brief This function computes a (quasi-random) list of bits to be taken
 * from the full descriptor. To speed the extraction, the function creates
 * a list of the samples that are involved in generating at least a bit (sampleList)
 * and a list of the comparisons between those samples (comparisons)
 * @param sampleList
 * @param comparisons The matrix with the binary comparisonsw
 * @param nbits The number of bits of the descriptor
 * @param pattern_size The pattern size for the binary descriptor
 * @param nchannels Number of channels to consider in the descriptor (1-3)
 * @note The function keeps the 18 bits (3-channels by 6 comparisons) of the
 * coarser grid, since it provides the most robust estimations
 */
void generateDescriptorSubsample(Mat& sampleList, Mat& comparisons, int nbits,
                                 int pattern_size, int nchannels) {

  int ssz = 0;
  for (int i = 0; i < 3; i++) {
    int gz = (i + 2)*(i + 2);
    ssz += gz*(gz - 1) / 2;
  }
  ssz *= nchannels;

  CV_Assert(ssz == 162*nchannels);
  CV_Assert(nbits <= ssz && "Descriptor size can't be bigger than full descriptor (486 = 162*3 - 3 channels)");

  // Since the full descriptor is usually under 10k elements, we pick
  // the selection from the full matrix.  We take as many samples per
  // pick as the number of channels. For every pick, we
  // take the two samples involved and put them in the sampling list

  Mat_<int> fullM(ssz / nchannels, 5);
  for (int i = 0, c = 0; i < 3; i++) {
    int gdiv = i + 2; //grid divisions, per row
    int gsz = gdiv*gdiv;
    int psz = divUp(2*pattern_size, gdiv);

    for (int j = 0; j < gsz; j++) {
      for (int k = j + 1; k < gsz; k++, c++) {
        fullM(c, 0) = i;
        fullM(c, 1) = psz*(j % gdiv) - pattern_size;
        fullM(c, 2) = psz*(j / gdiv) - pattern_size;
        fullM(c, 3) = psz*(k % gdiv) - pattern_size;
        fullM(c, 4) = psz*(k / gdiv) - pattern_size;
      }
    }
  }

  RNG rng(1024);
  const int npicks = divUp(nbits, nchannels);
  Mat_<int> comps = Mat_<int>(nchannels * npicks, 2);
  comps = 1000;

  // Select some samples. A sample includes all channels
  int count = 0;
  Mat_<int> samples(29, 3);
  Mat_<int> fullcopy = fullM.clone();
  samples = -1;

  for (int i = 0; i < npicks; i++) {
    int k = rng(fullM.rows - i);
    if (i < 6) {
      // Force use of the coarser grid values and comparisons
      k = i;
    }

    bool n = true;

    for (int j = 0; j < count; j++) {
      if (samples(j, 0) == fullcopy(k, 0) && samples(j, 1) == fullcopy(k, 1) && samples(j, 2) == fullcopy(k, 2)) {
        n = false;
        comps(i*nchannels, 0) = nchannels*j;
        comps(i*nchannels + 1, 0) = nchannels*j + 1;
        comps(i*nchannels + 2, 0) = nchannels*j + 2;
        break;
      }
    }

    if (n) {
      samples(count, 0) = fullcopy(k, 0);
      samples(count, 1) = fullcopy(k, 1);
      samples(count, 2) = fullcopy(k, 2);
      comps(i*nchannels, 0) = nchannels*count;
      comps(i*nchannels + 1, 0) = nchannels*count + 1;
      comps(i*nchannels + 2, 0) = nchannels*count + 2;
      count++;
    }

    n = true;
    for (int j = 0; j < count; j++) {
      if (samples(j, 0) == fullcopy(k, 0) && samples(j, 1) == fullcopy(k, 3) && samples(j, 2) == fullcopy(k, 4)) {
        n = false;
        comps(i*nchannels, 1) = nchannels*j;
        comps(i*nchannels + 1, 1) = nchannels*j + 1;
        comps(i*nchannels + 2, 1) = nchannels*j + 2;
        break;
      }
    }

    if (n) {
      samples(count, 0) = fullcopy(k, 0);
      samples(count, 1) = fullcopy(k, 3);
      samples(count, 2) = fullcopy(k, 4);
      comps(i*nchannels, 1) = nchannels*count;
      comps(i*nchannels + 1, 1) = nchannels*count + 1;
      comps(i*nchannels + 2, 1) = nchannels*count + 2;
      count++;
    }

    Mat tmp = fullcopy.row(k);
    fullcopy.row(fullcopy.rows - i - 1).copyTo(tmp);
  }

  sampleList = samples.rowRange(0, count).clone();
  comparisons = comps.rowRange(0, nbits).clone();
}

}
