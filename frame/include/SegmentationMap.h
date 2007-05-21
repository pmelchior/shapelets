#ifndef SEGMENTATION_MAP
#define SEGMENTATION_MAP

#include <Image.h>
#include <list>

/// Segmentation map class.
/// Methods for generating or working on segmentation maps. Most important:
/// - linkPixelsSetMap(): build list of connected pixels.
/// - findObjectPixels(): returns all pixels of a selected object.
/// - findHalo(): finds the halo arround an object typically unrecognized by 
///  a threshold detection.
/// 
/// \todo Grid maybe not necessary here?

class SegmentationMap : public Image<int> {
 public:
  /// Constructor for generating a new segmentation map from a FitsImage.
  SegmentationMap(const Image<double>& image);
  /// Constructor from the existing segmentation FITS file <tt>segMapFile</tt>.
  /// The <tt>image</tt> has to be defined on the same Grid as the given  
  /// segmentation map.
  SegmentationMap(std::string segMapFile, const Image<double>& image);
  /// Copy constructor
  SegmentationMap(const SegmentationMap& segMap);
  /// Copy operator
  void operator=(const SegmentationMap& segMap);
  /// Get the number of different objects present in the map.
  /// This is a brute-force method because it traverses the whole map.
  unsigned int getNumberOfObjects();
  /// Get a list of connected pixels which have values above (positive) or below 
  /// (negative) <tt>threshold</tt>.
  /// This is a Friend-of-Friend algorithm with a linking length of 1 pixel.\n
  /// It starts by putting <tt>startpixel</tt> into the <tt>pixellist</tt> and setting
  /// <tt>segMap(startpixel) = tag</tt> if <tt>imageData(startpixel) > (<) threshold</tt>.
  /// Then it performs these two steps for all neighbors of pixels in the list until no
  /// new pixels are found.
  void linkPixelsSetMap(std::list<unsigned int>& pixellist, unsigned int startpixel, int tag, double threshold, bool positive);
  /// Find the positive halo around the object which is typically unrecognized by a 
  /// threshold detection.
  /// Starting from <tt>corelist</tt> (already found pixels of object),
  /// this searches for positive (negative) pixel groups larger than 
  /// <tt>0.5*correlationLength^2</tt> (pixels) in an area 50% larger than 
  /// <tt>xmin..xmax, ymin..ymax</tt>. <tt>correlationLength</tt> can be obtained from
  /// parameters indicating seeing conditions or PSF of the instrument.\n
  /// For each of the positive pixels it determines the minimum distance to the object 
  /// and fills the distance into a histogram.  Since the halo should create a high number 
  /// of pixels with short distances, we define the halo as
  /// those pixels with smaller distances than the distance of the minimum after the first 
  /// peak in the histogram (see example below).
  /// \image html histogram_halo_small.png
  /// The coordinates <tt>xmin..xmax, ymin..ymax</tt> are updated to include the halo, 
  /// the halo pixels are set to <tt>objectnr</tt> in the segmentation map and added to <tt>corelist</tt>.\n
  /// For the found pixel groups not related to the halo, segmentation map is set to <tt>-1</tt>
  /// (positive) or <tt>-2</tt> (negative).
  void findHalo(unsigned int objectnr, std::list<unsigned int>& corelist, int& xmin, int& xmax, int& ymin, int& ymax, double correlationLength, double bg_mean, double bg_rms);
  /// Find list of pixels due to given object from segmentation map.
  /// It reads the segmentation map in the range of <tt>xmin..xmax,ymin..ymax</tt> only.
  void findObjectPixels(std::list<unsigned int>& pixellist, unsigned int objectnr, int xmin, int xmax, int ymin, int ymax);
  /// Removes all inner pixels of this object.
  /// From the <tt>pixellist</tt> (which should contain the list of pixels of the given object,
  /// filled by e.g. findObjectPixels()), all pixels are removed which have only neighbors 
  /// from the same object.
  void removeObjectInnerPixels(std::list<unsigned int>&  pixellist, unsigned int objectnr);
  /// Set all pixels on the border of the segment to <tt>tag</tt>.
  /// It "paints" a open rectangle of with a value of <tt>tag</tt> along the border 
  /// <tt>xmin..xmax,ymin..ymax</tt>.\n
  /// <b>Attention:</b> Since background pixel are considered to have <tt>segMap(pixel) == 0</tt>,
  /// changing values on the border effectively removes those pixels from the set of
  /// background pixels whenever <tt>tag!=0</tt>.
  void setSegmentBorder(int tag, int xmin, int xmax, int ymin, int ymax);
  /// Save segmentation map as FITS file.
  /// The extension has to given in the standard cfitsio way as <tt>filename[extension]</tt>.
  void save(std::string filename, std::map<std::string,std::string> keywords);
 private:
  NumVector<double>& data;
  void addFrameBorder(double factor, int& xmin, int& xmax, int& ymin, int& ymax);
  double distanceFromRim(std::list<unsigned int>& objectlist, unsigned int pixel);
  void cleanSegMapArea(int xmin, int xmax, int ymin, int ymax);
};

#endif
