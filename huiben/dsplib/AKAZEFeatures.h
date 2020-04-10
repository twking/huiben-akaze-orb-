/**
 * @file AKAZE.h
 * @brief Main class for detecting and computing binary descriptors in an
 * accelerated nonlinear scale space
 * @date Mar 27, 2013
 * @author Pablo F. Alcantarilla, Jesus Nuevo
 */

#ifndef __AKAZE_FEATURES_H__
#define __AKAZE_FEATURES_H__

/* ************************************************************************* */
// Includes
#include "AKAZEConfig.h"
#include "precomp.hpp"

namespace libakaze
{
using namespace libakaze;
/// A-KAZE nonlinear diffusion filtering evolution
template <typename MatType>
struct Evolution
{
  Evolution() {
    etime = 0.0f;
    esigma = 0.0f;
    octave = 0;
    sublevel = 0;
    sigma_size = 0;
    octave_ratio = 0.0f;
    border = 0;
  }

  template <typename T>
  explicit Evolution(const Evolution<T> &other) {
    size = other.size;
    etime = other.etime;
    esigma = other.esigma;
    octave = other.octave;
    sublevel = other.sublevel;
    sigma_size = other.sigma_size;
    octave_ratio = other.octave_ratio;
    border = other.border;

    other.Lx.copyTo(Lx);
    other.Ly.copyTo(Ly);
    other.Lt.copyTo(Lt);
    other.Lsmooth.copyTo(Lsmooth);
    other.Ldet.copyTo(Ldet);
  }

  MatType Lx, Ly;           ///< First order spatial derivatives
  MatType Lt;               ///< Evolution image
  MatType Lsmooth;          ///< Smoothed image, used only for computing determinant, released afterwards
  MatType Ldet;             ///< Detector response

  Size size;                ///< Size of the layer
  float etime;              ///< Evolution time
  float esigma;             ///< Evolution sigma. For linear diffusion t = sigma^2 / 2
  int octave;               ///< Image octave
  int sublevel;             ///< Image sublevel in each octave
  int sigma_size;           ///< Integer esigma. For computing the feature detector responses
  float octave_ratio;       ///< Scaling ratio of this octave. ratio = 2^octave
  int border;               ///< Width of border where descriptors cannot be computed
};

typedef Evolution<Mat> MEvolution;
typedef Evolution<Mat> UEvolution;//UMat
typedef std::vector<MEvolution> Pyramid;
typedef std::vector<UEvolution> UMatPyramid;

/* ************************************************************************* */
// AKAZE Class Declaration
class AKAZEFeatures {

private:

  AKAZEOptions options_;                ///< Configuration options for AKAZE
  Pyramid evolution_;        ///< Vector of nonlinear diffusion evolution

  /// FED parameters
  int ncycles_;                  ///< Number of cycles
  bool reordering_;              ///< Flag for reordering time steps
  std::vector<std::vector<float > > tsteps_;  ///< Vector of FED dynamic time steps
  std::vector<int> nsteps_;      ///< Vector of number of steps per cycle

  /// Matrices for the M-LDB descriptor computation
  libakaze::Mat descriptorSamples_;  // List of positions in the grids to sample LDB bits from.
  libakaze::Mat descriptorBits_;
  libakaze::Mat bitMask_;

  /// Scale Space methods
  void Allocate_Memory_Evolution();
  void Find_Scale_Space_Extrema(std::vector<Mat>& keypoints_by_layers);
  void Do_Subpixel_Refinement(std::vector<Mat>& keypoints_by_layers,std::vector<KeyPoint>& kpts);

  /// Feature description methods
  void Compute_Keypoints_Orientation(std::vector<KeyPoint>& kpts) const;

public:
  /// Constructor with input arguments
  AKAZEFeatures(const AKAZEOptions& options);
  void Create_Nonlinear_Scale_Space(InputArray img);
  void Feature_Detection(std::vector<KeyPoint>& kpts);
  void Compute_Descriptors(std::vector<KeyPoint>& kpts, OutputArray desc);
};

/* ************************************************************************* */
/// Inline functions
void generateDescriptorSubsample(libakaze::Mat& sampleList, libakaze::Mat& comparisons,
                                 int nbits, int pattern_size, int nchannels);

}

#endif
