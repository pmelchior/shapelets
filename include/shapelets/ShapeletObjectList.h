#ifndef SHAPELETOBJECTLIST_H
#define SHAPELETOBJECTLIST_H

#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>
#include <Typedef.h>
#include <shapelets/ShapeletObject.h>

/// Class for managing an ensemble of ShapeletObjects.
/// This class is meant for loading a set of ShapeletObject entities by giving a list
/// of SIFFile names as an ASCII file.\n
/// It is derived from <tt>std::vector<boost::shared_ptr<ShapeletObject> ></tt>, 
/// so effectively one has a vector of pointers to ShapeletObjects with the additional layer of the
/// <tt>boost::shared_ptr</tt> which ensures that created object are deleted when the last 
/// <tt>shared_ptr</tt> pointing to it is destroyed or reset.\n\n
/// For more details on <tt>std::vector</tt> and <tt>boost::shared_ptr</tt>, see
/// http://www.cppreference.com/cpplist/index.html and http://www.boost.org/libs/smart_ptr/shared_ptr.htm.\n\n
///
/// The elements of the vector can by accessed via an index or an iterator:
/// \code
/// // file.lst is generated by something like
/// // ls *.sif > file.lst
/// ShapeletObjectList sl("file.lst");
/// ShapletObjectList::iterator iter;
/// for (iter = sl.begin(); iter != sl.end(); iter++) {
///   ShapeletObject* sobj = (*iter);
///   ...
/// }
/// \endcode
///
/// Apart from that, one can select an arbitrary criterium for any ShapeletObject 
/// to select a subset of these objects by giving an appropriate
/// <tt>selectionFunction</tt>:
/// \code
/// bool selectionFunction(ShapeletObject& sobj) {
///   if (sobj.getDecompositionChiSquare() < 1) 
///     return 1;
///   else
///     return 0;
/// }
/// ShapeletObjectList sl("file.lst",&selectionFunction);
/// \endcode
/// In this example, all SIFFile names stored in <tt>file.lst</tt> are opened;
/// if their decomposition \f$\chi^2\f$ is less than 1, they will be included in
/// the list.


class ShapeletObjectList : public std::vector<boost::shared_ptr<ShapeletObject> > {
 public:
  /// Default contructor.
  ShapeletObjectList();
  /// Constructor for loading all files in <tt>listfile</tt>.
  ShapeletObjectList(std::string listfile);
  /// Constructor for loading a subset of files from <tt>listfile</tt>. 
  /// The ShapeletObject <tt>sobj</tt> entities have to fulfill the criterium 
  /// <tt>selectionFunction (sobj) == 1</tt> to be included in the list.
  ShapeletObjectList(std::string listfile, bool (* selectionFunction) (ShapeletObject&));
  /// Select a subset of ShapeletObject entities.
  /// The usage of <tt>selectionFunction</tt> is identical as in the
  /// constructor.
  ShapeletObjectList select(bool (* selectionFunction) (ShapeletObject&));
  /// Compute coefficient and scale size average of all SIFs in the list.
  /// <tt>std_mean</tt> is the standard deviation of the mean of the coefficient set.
  void average(NumMatrix<data_t>& mean, NumMatrix<data_t>& std_mean, data_t& beta);
  /// Compute coefficient and scale size average of all ShapeletObject entities in the list.
  /// <tt>std_mean</tt> is the standard deviation of the mean of the coefficient set.
  ///
  /// One can weight the ShapeletObject entities according to an arbitrary <tt>weightFunction</tt>. For example:
  /// \code
  /// data_t weight(ShapeletObject& sobj) {
  ///   return sobj.getShapeletFlux();
  /// }
  /// NumMatrix<data_t> mean, std_mean;
  /// data_t beta;
  /// sl.average(mean,std_mean,beta,&weight);
  /// \endcode
  ///
  /// This returns the (flux-) weighted mean
  /// \f[\langle x\rangle = \frac{ \sum_{i=1}^N w_i x_i}{\sum_{i=1}^N w_i}\f]
  /// and its standard deviation 
  /// \f[s_{\langle x\rangle} = \frac{ \sum_{i=1}^N{w_i} }{(\sum_{i=1}^N{w_i})^2 - \sum_{i=1}^N{{w_i}^2}} \sum_{i=1}^N{{w_i}(x_i - \langle x\rangle)^2}\f]
  /// which is computed as 
  /// \f[s_{\langle x\rangle} = \frac{\sum_{i=1}^N{w_i {x_i}^2} \sum_{i=1}^N{w_i} - (\sum_{i=1}^N{w_i x_i})^2} {(\sum_{i=1}^N{w_i})^2 - \sum_{i=1}^N{{w_i}^2}}\f]
  /// where \f$N\f$ is <tt>sl.size()</tt> and \f$x\f$ is any available coefficient in <tt>sl</tt>.
  void average(NumMatrix<data_t>& mean, NumMatrix<data_t>& std_mean, data_t& beta, data_t (* weightFunction) (ShapeletObject&));

 private:
  void readListFile(std::string listfile, bool (* selectionFunction) (ShapeletObject&));
  bool checkSIFFile(ShapeletObject&, std::string );
};
#endif
