#ifndef SHAPELENS_H
#define SHAPELENS_H

#include "Typedef.h"
#include "ShapeLensConfig.h"

#include "frame/Frame.h"
#include "frame/SExFrame.h"

#include "shapelets/ShapeletObject.h"
#include "shapelets/SIFFile.h"
#include "shapelets/ShapeletObjectList.h"

#include "modelfit/SourceModel.h"

/// Namespace for ShapeLens++
namespace shapelens {

/** \mainpage ShapeLens++ Documentation
\section intro Introduction
The Shapelets formalism is a recent approach to image analysis, especially suited for
observational astronomy.\n
In principle one decomposes an arbitrary 2D function into a series
of localized orthogonal basis functions, which are called 'shapelets'. The specialty of this
basis set is that it consists of weighted Hermite polynomials, which correspond to
perturbations arround a Gaussian. The 'shapelets' are also the basis function of the 2D
quantum harmonic oscillator, and thus allow the use of the operator formalism developed
for this problem: transformations like rotations, translations etc. can be described
by combinations of lowering and raising operators, acting on the 'shapelets'.\n
Since galaxy objects or their constituents have a localized appearance, they are conveniently
described by using a basis set of localized functions. Thus the number of coefficients necessary
to describe a galaxy object well will generally be low. This makes the Shapelets formalism
interesting for data reduction and storage in the first place. Furthermore, as it is possible to
to associate physical quantities of the galaxies to functions of their shapelet coefficients,
one can do the physical analysis in the much smaller shapelet space instead of the real
space, thus saving memory and computation time.

\section refs References
- Refregier A., MNRAS, 2003, 338, 35 (later called: Paper I)
- Refregier A. & Bacon D., MNRAS, 2003, 338, 48 (Paper II)
- Massey R. & Refregier A.,MNRAS, 2005, 363, 197 (Paper III)
- Melchior, P. et al., A&A, 2007, 463, 1215 (Paper IV)

\section library Libary Overview
Although not every line of code is pure C++, the approach is purely Object-oriented. All 
necessary methods for constructing an analysis pipeline are public, well-documented 
and named in a sensible fashion.
This renders the development of a specialzed version of this code straightforward.\n\n
The workflow is the following:
- Reading image information, segmenting images into individual object and storing all
necessary information of the pixelized object in one structure (an instance of the 
Object class) is the purpose of the Frame class. 
How this can be generalized to read e.g. SExtractor outputs
is demonstrated by the SExFrame class.
- The Object is the passed to the ShapeletObject class, which does the actual shapelet
decomposition by calling OptimalDecomposite and Decomposite. The decomposition
result can then be queried and save to a FITS-compatible SIFFile.
- All further analysis is then done on a single ShapeletObject or on an ensemble of them,
organized as a ShapeletObjectList.

\section externals Additional Libraries
- cfitsio: library to work with FITS files, in C (http://heasarc.gsfc.nasa.gov/docs/software/fitsio/fitsio.html).
- GNU Scientific Library (GSL): Free and well tested implementation of common
(and also not so common) mathematical functions and procedures in C
(http://www.gnu.org/software/gsl/).
- boost: various extensions to standard C++ (http://www.boost.org/)
- numla: library for doing vector and matrix operations conveniently and fast (https://www.ita.uni-heidelberg.de/internal/projects/numla)

\author Peter Melchior (pmelchior at ita dot uni-heidelberg dot de)
*/

} // end namespace
#endif
