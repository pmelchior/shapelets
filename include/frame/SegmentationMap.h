#ifndef SEGMENTATION_MAP
#define SEGMENTATION_MAP

#include <set>
#include <Typedef.h>
#include <frame/Image.h>

/// Segmentation map class.
/// Methods for working on segmentation maps.

class SegmentationMap : public Image<long> {
 public:
  /// Defautl constructor.
  SegmentationMap();
  /// Constructor from the existing segmentation FITS file <tt>segMapFile</tt>.
  /// The <tt>image</tt> has to be defined on the same Grid as the given  
  /// segmentation map.
  SegmentationMap(std::string segMapFile);
  /// Get the number of different objects present in the map.
  /// This is a brute-force method because it traverses the whole map.
  unsigned long getNumberOfObjects();
  /// Find set of pixels due to given object from segmentation map.
  /// It reads the segmentation map in the range of <tt>xmin..xmax,ymin..ymax</tt> only.
  void findObjectPixels(std::set<unsigned long>& pixellist, unsigned long objectnr, long xmin, long xmax, long ymin, long ymax);
  /// Set all pixels on the border of the segment to <tt>tag</tt>.
  /// It "paints" a open rectangle of with a value of <tt>tag</tt> along the border 
  /// <tt>xmin..xmax,ymin..ymax</tt>.\n
  /// <b>Attention:</b> Since background pixel are considered to have <tt>segMap(pixel) == 0</tt>,
  /// changing values on the border effectively removes those pixels from the set of
  /// background pixels whenever <tt>tag!=0</tt>.
  void setSegmentBorder(int tag, long xmin, long xmax, long ymin, long ymax);
  /// Clean the segmentation map within the given region from already detected objects.
  void cleanSegMapArea(long xmin, long xmax, long ymin, long ymax);
};

#endif
