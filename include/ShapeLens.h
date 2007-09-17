#ifndef SHAPELENS_H
#define SHAPELENS_H

/** \mainpage ShapeLens++ Documentation
\section Introduction
The Shapelets formalism is a recent approach to image analysis, especially suited for
observational astronomy.\n
In principle one decomposes an arbitrary 2D function into a series
of localized orthogonal basis functions, which are called 'shapelets'. The specialty of this
basis set is that it consists of weighted Hermite polynomials, which correspond to
perturbations arround a Gaussian. The 'shapelets' are also the basis function of the 2D
quantum harmonic oscillator, and thus allow the use of the operator formalism developed
for this problem. E.g. transformations as rotations, translations etc. can be described
by combinations of lowering and raising operators, acting on the 'shapelets'.\n
Since galaxy objects or their constituents have a localized appearance, they are conveniently
described by using a basis set of localized functions. Thus the number of coefficients needed
to describe a galaxy object well will generally be low. This makes the Shapelets formalism
interesting for data reduction and storage in the first place. Further, as it is possible to
to associate physical quantities of the galaxies to functions of their shapelet coefficients,
one can do the physical analysis in the much smaller shapelet space instead of the real
space, thus saving memory and computation time.
\section References
- Refregier A., MNRAS, 2003, 338, 35 (later called: Paper I)
- Refregier A. & Bacon D., MNRAS, 2003, 338, 48 (Paper II)
- Massey R. & Refregier A.,MNRAS, 2005, 363, 197 (Paper III)
- Melchior, P. et al., A&A, 2007, 463, 1215 (Paper IV)

\section Additional Libraries
- cfitsio: library to work with FITS files, in C (http://heasarc.gsfc.nasa.gov/docs/software/fitsio/fitsio.html).
- GNU Scientific Library (GSL): Free and well tested implementation of common
(and also not so common) mathematical functions and procedures in C
(http://www.gnu.org/software/gsl/).
- boost: various extensions to standard C++ (http://www.boost.org/)
- numla: library for doing vector and matrix operations conveniently and fast (https://www.ita.uni-heidelberg.de/internal/projects/numla)


\author Peter Melchior (pmelchior at ita dot uni-heidelberg dot de)
*/

#include <IO.h>
#include <ShapeLensConfig.h>
#include <History.h>

#include <frame/Frame.h>
#include <frame/SExFrame.h>

#include <shapelets/ShapeletObject.h>
#include <shapelets/ShapeletObjectList.h>

#include <lensing/LensingEstimator.h>

#endif
