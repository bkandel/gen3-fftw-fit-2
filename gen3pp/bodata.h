#ifndef _BODATA_H_
#define _BODATA_H_
//
// bodata.h
// Provides an interface for multi-dimensional data (vector, matrix, 3d-matrix, 4d-matrix etc)
// Author: SauravP
// 
#include <vector>

//
// The main data class.
// * This class can be used for many dimensional data
// * The data is encoded as:
//   * index = i + axes[0]*(j + axes[1]*(k + ....
template<typename T>
struct BOData {
    std::vector<T> data;
    std::vector<long> axes;
};

#endif
