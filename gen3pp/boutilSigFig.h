#ifndef _BOUTIL_H_
#define _BOUTIL_H_
//
// boutil.h
// Provides utilities for multi-dimensional data (vector, matrix, 3d-matrix, 4d-matrix etc)
// Author: SauravP
// 
#include <valarray>
#include <vector>
#include <utility>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <gen3pp/bodata.h>

//
// Print some generic stats
//
template<typename T>
int print(BOData<T>& bd) {
    std::cout << "BOData stats:" << std::endl;
    std::cout << "  size: " << bd.data.size() << std::endl;
    for (unsigned int i=0; i<bd.axes.size(); ++i) {
        std::cout << "  axes: " << bd.axes[i] << std::endl;
    }
    std::cout << "  Max: " << *max_element(bd.data.begin(), bd.data.end()) << std::endl; 
    std::cout << "  Min: " << *min_element(bd.data.begin(), bd.data.end()) << std::endl; 
    return 0;
}

template<typename T>
int resetBOData2(BOData<T>& d, unsigned int isize, unsigned int jsize, T x) {
    d.data.clear();
    d.axes.clear();
    d.data.resize(isize*jsize, x);
    d.axes.push_back(isize);
    d.axes.push_back(jsize);
    return 0;
}

template<typename T>
int allocateBOData(const std::vector<long int>& axes, BOData<T>& data) {
    data.data.clear();
    data.axes.clear();
    unsigned int dim = axes.size();
    long n = 1;
    
    data.axes.resize(dim);
    for (unsigned int i=0; i<dim; ++i) {
        data.axes[i] = axes[i];
        n *= data.axes[i];
    }
    data.data.resize(n);
    return 0;
}     
    

//
// Add BOData x to to
//
template<typename T>
int addto(const BOData<T>& x, BOData<T>& t) {
    for (int i=0; i<x.axes.size(); ++i) {
        if (x.axes[i] != t.axes[i]) {
            std::cerr << "Dimensions not equal" << std::endl;
            return 1;
        }
    } 
    
    for (int i=0; i<t.data.size(); ++i) {
        t.data[i] += x.data[i];
    }
}


//
// Given 3D data "from" substract 2D data "by" from each of its frames
// 
template<typename T> 
int subtract(const BOData<T>& by, BOData<T>& from) {
    if (by.axes[0] != from.axes[0] || by.axes[1] != from.axes[1]) {
        std::cerr << "Dimensions not equal" << std::endl;
        return 1;
    } 
    if (by.axes.size() == 2 && from.axes.size() == 2) {
        for (unsigned int i=0; i<from.data.size(); ++i) {
            from.data[i] -= by.data[i];
        }
        return 0;
    }

    {
        unsigned int isize = from.axes[0];
        unsigned int jsize = from.axes[1];
        unsigned int ksize = from.axes[2];
        for (unsigned int i=0; i<isize; ++i) {
            for (unsigned int j=0; j<jsize; ++j) {
                T by_t = by.data[i+isize*j];
                for (unsigned int k=0; k<ksize; ++k) {
                    unsigned int idx  = i + isize*(j + jsize*k);
                    //from.data[idx] -= by.data[nidx];
                    from.data[idx] -= by_t;
                }
            }
        }
        return 0;
    }
}

//
// Swap two elements in the data
//
template<typename T>
int swap(BOData<T>& d, unsigned int i, unsigned int j) {
    T temp = d.data[i];
    d.data[i] = d.data[j];
    d.data[j] = temp;
} 

//
//  Transpose of a two dimensional matrix.
//
template<typename T>
int transpose(BOData<T>& d) {
    //
    // do nothing if it is not a two dimensional matrix
    if (d.axes.size() != 2) {
        return 0;
    }

    for (unsigned int i=0; i<d.axes[0]; ++i) {
        for (unsigned int j=i+1; j<d.axes[1]; ++j) {
            unsigned int ii = i + d.axes[0]*j;
            unsigned int jj = j + d.axes[1]*i;
            swap(d, ii, jj);
        }
    }
    unsigned int i = d.axes[0];
    d.axes[0] = d.axes[1];
    d.axes[1] = i;
    return 0;
}

//
// Flip about a horizontal axis
//
template<typename T>
int flipUD(BOData<T>& d) {
    //
    // do nothing if it is not a two dimensional matrix
    if (d.axes.size() != 2) {
        return 0;
    }

    for (unsigned int i=0; i<d.axes[0]; ++i) {
        for (unsigned int j=0; j<d.axes[1]/2; ++j) {
            unsigned int ii = i + d.axes[0]*j;
            unsigned int jj = i + d.axes[0]*(d.axes[1]-1-j);
            swap(d, ii, jj);
        }
    }
    
    return 0;
} 
       
//
// Flip about a vertical axis
//            
template<typename T>
int flipLR(BOData<T>& d) {
    //
    // do nothing if it is not a two dimensional matrix
    if (d.axes.size() != 2) {
        return 0;
    }

    for (unsigned int i=0; i<d.axes[0]/2; ++i) {
        for (unsigned int j=0; j<d.axes[1]; ++j) {
            unsigned int ii = i + d.axes[0]*j;
            unsigned int jj = (d.axes[0]-1-i) + d.axes[0]*j;
            swap(d, ii, jj);
        }
    }
    
    return 0;
} 

//
// Rotate anti-clockwise
//
template<typename T> 
int rotateAntiClockwise(BOData<T>& d) {
    transpose(d);
    flipUD(d);
} 

template<typename T> 
int rotate90(BOData<T>& d) {
    transpose(d);
    flipUD(d);
} 
    
    
//
// Rotate clockwise
//
template<typename T> 
int rotateClockwise(BOData<T>& d) {
    transpose(d);
    flipLR(d);
} 

//
// Get the mean along a given dimension
//
template<typename T>
int getMean(BOData<T>& data, unsigned int c, BOData<T>& mean) {
    long isize = data.axes[0];
    long jsize = data.axes[1];
    long ksize = data.axes[2];
    unsigned int k, j, i;
    unsigned int *i1, *i2, *ix;
    unsigned int stride;

    mean.data.clear();
    mean.axes.clear();
    switch (c) {
        case 1:
            ix = &i;
            i1 = &j;
            i2 = &k;
            stride = jsize;
            mean.data.resize(jsize*ksize, 0.0);
            mean.axes.push_back(jsize);
            mean.axes.push_back(ksize);
            break;
        case 2:
            i1 = &i;
            ix = &j;
            i2 = &k;
            stride = isize;
            mean.data.resize(isize*ksize, 0.0);
            mean.axes.push_back(isize);
            mean.axes.push_back(ksize);
            break;
        case 3:
            i1 = &i;
            i2 = &j;
            ix = &k;
            stride = isize;
            mean.data.resize(isize*jsize, 0.0);
            mean.axes.push_back(isize);
            mean.axes.push_back(jsize);
            break;
    }
    for (k=0; k<ksize; ++k) {
        for (j=0; j<jsize; ++j) {
            for (i=0; i<isize; ++i) {
                unsigned int idx  = i + isize*(j + jsize*k);
                unsigned int nidx = *i1 + *i2*stride;
                //mean.data[nidx] += (static_cast<double>(data.data[idx]) - mean.data[nidx])/(*ix+1);
                mean.data[nidx] += (data.data[idx] - mean.data[nidx])/(*ix+1);
            }
        }
    }
    return 0;
}

template<typename T>
int getVar(BOData<T>& data, unsigned int c, BOData<T>& var) {

    BOData<T> mean;
    getMean(data, c, mean);
    
    long isize = data.axes[0];
    long jsize = data.axes[1];
    long ksize = data.axes[2];
    unsigned int k, j, i;
    unsigned int *i1, *i2, *ix;
    unsigned int stride;

    var.data.clear();
    var.axes.clear();
    switch (c) {
        case 1:
            ix = &i;
            i1 = &j;
            i2 = &k;
            stride = jsize;
            var.data.resize(jsize*ksize, 0.0);
            var.axes.push_back(jsize);
            var.axes.push_back(ksize);
            break;
        case 2:
            i1 = &i;
            ix = &j;
            i2 = &k;
            stride = isize;
            var.data.resize(isize*ksize, 0.0);
            var.axes.push_back(isize);
            var.axes.push_back(ksize);
            break;
        case 3:
            i1 = &i;
            i2 = &j;
            ix = &k;
            stride = isize;
            var.data.resize(isize*jsize, 0.0);
            var.axes.push_back(isize);
            var.axes.push_back(jsize);
            break;
    }

    for (k=0; k<ksize; ++k) {
        for (j=0; j<jsize; ++j) {
            for (i=0; i<isize; ++i) {
                unsigned int idx  = i + isize*(j + jsize*k);
                unsigned int nidx = *i1 + *i2*stride;
                const double delta = data.data[idx] - mean.data[nidx];
                var.data[nidx] += (delta*delta - var.data[nidx])/(*ix+1);
            }
        }
    }
    return 0;
}

template<typename T>
int getSD(BOData<T>& data, unsigned int c, BOData<T>& sd) {
    getVar(data, c, sd);
    for (unsigned int i=0; i<sd.data.size(); ++i) {
        sd.data[i] = sqrt(sd.data[i]);
    }
    return 0;
}
 
template<typename T>
int smooth2D(const BOData<T>& data, unsigned int c, BOData<T>& smooth) {
    if (data.axes.size() != 2) {
        std::cerr << "Wrong dimensions of data in smooth2D(): " << data.axes.size() << std::endl;
        return 1;
    }

    if (((c%2)!=1) && (c<=1)) {
        std::cerr << "Smoothing size is not an odd number greater than 1: " << c << std::endl;
        return 2;
    }

    unsigned int dc = (c-1)/2;

    smooth.data.clear();
    smooth.axes.clear();
    smooth.data.resize(data.data.size(), 0);
    smooth.axes = data.axes;

    for (unsigned int i=0; i<data.axes[0]; ++i) {
        for (unsigned int j=0; j<data.axes[1]; ++j) {
            unsigned int xmin = (i <= dc) ? 0 : i-dc;
            unsigned int ymin = (j <= dc) ? 0 : j-dc;
            unsigned int xmax = (i+dc >= data.axes[0] -1) ? data.axes[0] : i+dc;
            unsigned int ymax = (j+dc >= data.axes[1] -1) ? data.axes[1] : j+dc;
            unsigned int n    = 0;
            for (unsigned int x=xmin; x<xmax; ++x) {
                for (unsigned int y=ymin; y<ymax; ++y) {
                    smooth.data[i+smooth.axes[0]*j] += 
                        (data.data[x + data.axes[0]*y] - smooth.data[i+smooth.axes[0]*j])/(++n);
                }
            }
        }
    }
    return 0;
}

template<typename T>
int write2file2(const std::string filename, BOData<T>& bd) {
    std::ofstream ofs;
    ofs.open(filename.c_str());
    
    unsigned int isize = bd.axes[0];
    unsigned int jsize = bd.axes[1];
    for (unsigned int j=0; j<jsize; ++j) {
        for (unsigned int i=0; i<isize; ++i) {
            ofs << bd.data[i+isize*j] << " ";
        }
    }
    ofs.close();
    return 0;
}

template<typename T>
int write2col(const std::string filename, const BOData<T>& bd) {
    std::ofstream ofs;
    ofs.open(filename.c_str());

    for (unsigned int i=0; i<bd.data.size(); ++i) {
        ofs << bd.data[i] << std::endl;
    }
    ofs.close();
    return 0;
}

template<typename T>
int write2gim(const std::string filename, const BOData<T>& bd) {
    std::ofstream ofs;
    ofs.open(filename.c_str());

    for (unsigned int i=0; i<bd.axes.size(); ++i) {
        ofs << bd.axes[i] << " ";
    }
    ofs << std::endl;
    for (unsigned int i=0; i<bd.data.size(); ++i) {
        ofs << bd.data[i] << " ";
    }
    ofs.close();
    return 0;
}

template<typename T>
int readgim(const std::string filename, BOData<T>& bd) {
    std::ifstream ifs;
    ifs.open(filename.c_str());
    std::string line;
    if (!ifs) {
        std::cerr << "Unable to read GIM file: " << filename << std::endl;
        return 1;
    }
    unsigned int size = 1;
    {
        getline(ifs, line);
        bd.axes.clear();
        std::stringstream ss(line); 
        unsigned int dummy;
        while (ss >> dummy) {
            bd.axes.push_back(dummy);
            size *= dummy;
            //std::cerr << dummy << " ";
        }
    }
    //std::cerr << std::endl << "size: " << size << std::endl;
    T dummy;
    bd.data.resize(size);
    for (unsigned int i=0; i<size; ++i) {
        ifs >> bd.data[i];
        //bd.data[i] = dummy;
    }
    return 0;
} 

//
// Write to BOD format
//
template<typename T>
int writeBOD(const std::string filename, BOData<T>& bd) {
    std::ofstream ofs;
    ofs.open(filename.c_str());

    for (unsigned int i=0; i<bd.axes.size(); ++i) {
        ofs << bd.axes[i] << " ";
    }
    ofs << std::endl;
    
    for (unsigned int j=0; j<bd.data.size(); ++j) {
        ofs << bd.data[j] << std::endl;
    }
    ofs.close();
    return 0;
}

//
// Read BOD format
//
template<typename T>
int readBOD(const std::string filename, BOData<T>& bd) {
    std::ifstream ifs;
    ifs.open(filename.c_str());

    std::string line;
    getline(ifs, line);
    std::stringstream ss(line);

    bd.data.clear();
    bd.axes.clear();
    unsigned int temp;
    unsigned int size = 1;
    while (ss >> temp) {
        bd.axes.push_back(temp);
        size *= temp; 
    }
    bd.data.resize(size);
    T dat; 
    for (unsigned int i=0; i<size; ++i) {
        ifs >> dat;
        bd.data[i] = dat;
    }
    ifs.close();
    return 0;
}

template<typename T>
bool CmpDIndx(const std::pair<T, unsigned int>& a, const std::pair<T, unsigned int>&b) {
    return a.first < b.first;
}

template<typename T>
int centroid(const BOData<T>& data, unsigned int top, std::vector<std::pair<unsigned int, unsigned int> >& d2) {

    std::vector<std::pair<T, unsigned int> > vDIndx;
    for (unsigned int i=0; i<data.data.size(); ++i) {
        vDIndx.push_back(std::pair<T, unsigned int>(data.data[i], i));  
    }
    
    //std::sort(vDIndx.begin(), vDIndx.end(), CmpDIndx);
    std::sort(vDIndx.rbegin(), vDIndx.rend());

    unsigned int csize = data.axes[0];
    unsigned int rsize = data.axes[1];

    d2.clear();
    double rmean=0, cmean=0;
    for (unsigned int i=0; i<top; ++i) {
        unsigned int r = vDIndx[i].second / csize;
        unsigned int c = vDIndx[i].second % csize; 
        rmean += (static_cast<double>(r) - rmean)/(i+1);
        cmean += (static_cast<double>(c) - cmean)/(i+1);
        d2.push_back(std::pair<unsigned int, unsigned int>(static_cast<unsigned int>(rmean+0.5), static_cast<unsigned int>(cmean+0.5)));
    }

    return 0;
}

template<typename T>
int centroid(const BOData<T>& data, T topfraction, unsigned int& row_centroid, unsigned int& col_centroid) {

    if (data.axes.size() != 2) {
        std::cerr << "centroid: dim=" << data.axes.size() << ". This works only for 2D data" << std::endl;
        return 1;
    }

    // create an indexed vector of the data
    std::vector<std::pair<T, unsigned int> > vDIndx;
    for (unsigned int i=0; i<data.data.size(); ++i) {
        vDIndx.push_back(std::pair<T, unsigned int>(data.data[i], i));  
    }
    
    //std::sort(vDIndx.begin(), vDIndx.end(), CmpDIndx);
    //
    // sort it in decreasing order
    std::sort(vDIndx.rbegin(), vDIndx.rend());

    unsigned int csize = data.axes[0];
    unsigned int rsize = data.axes[1];
    unsigned int tsize = csize*rsize;
    unsigned int top   = static_cast<unsigned int>(static_cast<T>(tsize)*topfraction); // the topfraction elements

    T rmean=0, cmean=0;
    for (unsigned int i=0; i<top; ++i) {
        unsigned int r = vDIndx[i].second / csize;
        unsigned int c = vDIndx[i].second % csize; 
        rmean += (static_cast<T>(r) - rmean)/(i+1);
        cmean += (static_cast<T>(c) - cmean)/(i+1);
        //d2.push_back(std::pair<unsigned int, unsigned int>(static_cast<unsigned int>(rmean+0.5), static_cast<unsigned int>(cmean+0.5)));
    }
    row_centroid = static_cast<unsigned int>(rmean+0.5);
    col_centroid = static_cast<unsigned int>(cmean+0.5);

    return 0;
}


template<typename T>
int getMax(BOData<T>& data, BOData<T>& maxf) {
    long isize = data.axes[0];
    long jsize = data.axes[1];
    long ksize = data.axes[2];
    maxf.data.clear();
    maxf.axes.clear();

    maxf.axes.push_back(isize);
    maxf.axes.push_back(jsize);
    maxf.resize(isize*jsize, 0);

    for (unsigned int k=0; k<ksize; ++k) {
        for (unsigned int j=0; j<jsize; ++j) {
            for (unsigned int i=0; i<isize; ++i) {
                unsigned int idx  = i + isize*(j + jsize*k);
                maxf[i+isize*j] = (maxf[i+isize*j] < data.data[idx]) ? data.data[idx] : maxf.data[i+isize*j];
            }
        }
    }
    return 0;
}

template<typename T>
int getMin(BOData<T>& data, BOData<T>& minf) {
    long isize = data.axes[0];
    long jsize = data.axes[1];
    long ksize = data.axes[2];
    minf.data.clear();
    minf.axes.clear();

    minf.axes.push_back(isize);
    minf.axes.push_back(jsize);
    minf.resize(isize*jsize, 1e300);

    for (unsigned int k=0; k<ksize; ++k) {
        for (unsigned int j=0; j<jsize; ++j) {
            for (unsigned int i=0; i<isize; ++i) {
                unsigned int idx  = i + isize*(j + jsize*k);
                minf[i+isize*j] = (minf[i+isize*j] > data.data[idx]) ? data.data[idx] : minf.data[i+isize*j];
            }
        }
    }
    return 0;
}

template<typename T>
int subtractFrom(BOData<T>& by, BOData<T>& from) {

    if (by.axes.size() != from.axes.size())
        return 1;
    for (unsigned int i=0; i<by.axes.size(); ++i) {
        if (by.axes[i] != from.axes.axes[i])
            return 2;
    }

    for (unsigned int i=0; i<by.data.size(); ++i)
        from.data[i] -= by.data[i];

    return 0;
}

template<typename T>
int addTo(BOData<T>& by, BOData<T>& to) {

    if (by.axes.size() != to.axes.size())
        return 1;
    for (unsigned int i=0; i<by.axes.size(); ++i) {
        if (by.axes[i] != to.axes.axes[i])
            return 2;
    }

    for (unsigned int i=0; i<by.data.size(); ++i)
        to.data[i] += by.data[i];

    return 0;
}

template<typename T>
int getDiff(BOData<T>& a, BOData<T>& b, BOData<T> diff) {

    if (a.axes.size() != b.axes.size())
        return 1;
    for (unsigned int i=0; i<b.axes.size(); ++i) {
        if (a.axes[i] != b.axes.axes[i])
            return 2;
    }

    diff.data.clear();
    diff.axes.clear();
    diff.data.resize(a.data.size());
    for (unsigned int i=0; i<b.axes.size(); ++i) {
        diff.axes.push_back(b.axes[i]);
    }

    for (unsigned int i=0; i<b.data.size(); ++i)
        diff.data[i] = a.data[i]-b.data[i];

    return 0;
}

template<typename T>
int getSum(BOData<T>& a, BOData<T>& b, BOData<T> sum) {

    if (a.axes.size() != b.axes.size())
        return 1;
    for (unsigned int i=0; i<b.axes.size(); ++i) {
        if (a.axes[i] != b.axes.axes[i])
            return 2;
    }

    sum.data.clear();
    sum.axes.clear();
    sum.data.resize(a.data.size());
    for (unsigned int i=0; i<b.axes.size(); ++i) {
        sum.axes.push_back(b.axes[i]);
    }

    for (unsigned int i=0; i<b.data.size(); ++i)
        sum.data[i] = a.data[i]+b.data[i];

    return 0;
}

template<typename T>
int getMaxMinMean(BOData<T>& data, unsigned int c, BOData<T>& mmmean) {
    getMax(data, mmmean);
    BOData<T> minf;
    getMin(data, minf);

    getSubstract2(minf, mmmean);
    return 0;
}

template<typename T>
int divideBy(T div, BOData<T>& data) {
    for (unsigned int i=0; i<data.data.size(); ++i)
        data.data[i] /= div;
    return 0;
}

#endif
