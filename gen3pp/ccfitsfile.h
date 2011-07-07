#ifndef _CCFITSFILE_H_
#define _CCFITSFILE_H_
#include <string>
#include <valarray>
#include <vector>
#include <CCfits/CCfits>

/*******************************************************************************
	Class Name:
	Member names:
	Member Functions: 
	Base Classes:
	Derived Classes:
	Description
********************************************************************************/
//! The FITSRead class.
/*!
	This class is used for reading in the FITS files.
*/
class CCFITSfile {
	public:
        CCFITSfile(){}

        template<typename T>
        int read(const std::string& fileName, BOData<T>& bd) {
            CCfits::FITS pInfile(fileName.c_str(), CCfits::Read, true);
            CCfits::PHDU& image = pInfile.pHDU();
            image.readAllKeys();

            bd.data.clear();
            std::valarray<T> va;
            image.read(va);
            bd.data.resize(va.size());
            std::copy(&va[0], &va[0]+va.size(), bd.data.begin());
            va.resize(0);
            
            bd.axes.clear();
            for (unsigned int i=0; i<image.axes(); ++i) {
                bd.axes.push_back(image.axis(i));
            }
            return 0;
        }

};
#endif
