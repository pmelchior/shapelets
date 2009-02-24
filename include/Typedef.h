#ifndef TYPEDEF_H
#define TYPEDEF_H

/// Floating point format for (almost) all data used in ShapeLens++.
typedef double data_t;

/// Grid position type.
typedef int grid_t; 
template <class T> class GridT;
/// Grid type.
typedef GridT<grid_t> Grid;

#endif
