/**
 * @file AKAZEConfig.h
 * @brief AKAZE configuration file
 * @date Feb 23, 2014
 * @author Pablo F. Alcantarilla, Jesus Nuevo
 */

#ifndef __AKAZE_CONFIG_H__
#define __AKAZE_CONFIG_H__

namespace libakaze
{
/* ************************************************************************* */
/// AKAZE configuration options structure
struct AKAZEOptions {
	enum
	{
		DESCRIPTOR_KAZE_UPRIGHT = 2, ///< Upright descriptors, not invariant to rotation
		DESCRIPTOR_KAZE = 3,
		DESCRIPTOR_MLDB_UPRIGHT = 4, ///< Upright descriptors, not invariant to rotation
		DESCRIPTOR_MLDB = 5
	};
	enum
	{
		DIFF_PM_G1 = 0,
		DIFF_PM_G2 = 1,
		DIFF_WEICKERT = 2,
		DIFF_CHARBONNIER = 3
	};
    AKAZEOptions()
        : omax(1)
        , nsublevels(1)
        , img_width(0)
        , img_height(0)
        , soffset(1.6f)
        , derivative_factor(1.5f)
        , sderivatives(1.0)
        , diffusivity(DIFF_PM_G2)

        , dthreshold(0.0007f)
        , min_dthreshold(0.00001f)

        , descriptor(DESCRIPTOR_MLDB)
        , descriptor_size(256)
        , descriptor_channels(3)
        , descriptor_pattern_size(10)

        , kcontrast(0.001f)
        , kcontrast_percentile(0.7f)
        , kcontrast_nbins(300)
    {
    }

    int omax;                       ///< Maximum octave evolution of the image 2^sigma (coarsest scale sigma units)
    int nsublevels;                 ///< Default number of sublevels per scale level
    int img_width;                  ///< Width of the input image
    int img_height;                 ///< Height of the input image
    float soffset;                  ///< Base scale offset (sigma units)
    float derivative_factor;        ///< Factor for the multiscale derivatives
    float sderivatives;             ///< Smoothing factor for the derivatives
    int diffusivity;   ///< Diffusivity type

    float dthreshold;               ///< Detector response threshold to accept point
    float min_dthreshold;           ///< Minimum detector threshold to accept a point

    int descriptor;     ///< Type of descriptor
    int descriptor_size;            ///< Size of the descriptor in bits. 0->Full size
    int descriptor_channels;        ///< Number of channels in the descriptor (1, 2, 3)
    int descriptor_pattern_size;    ///< Actual patch size is 2*pattern_size*point.scale

    float kcontrast;                ///< The contrast factor parameter
    float kcontrast_percentile;     ///< Percentile level for the contrast factor
    int kcontrast_nbins;            ///< Number of bins for the contrast factor histogram
};

}

#endif
