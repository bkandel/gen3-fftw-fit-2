#ifndef _FFTWFIT_H_
#define _FFTWFIT_H_
#include <cmath>
#include <vector>
#include <iostream>
#include <fftw3.h>
#include <gen3pp/bodata.h>
#include <gen3pp/boutil.h>

//template<typename T>
//int resetBOData2(BOData<T>& d, unsigned int isize, unsigned int jsize, T x) {
//    d.data.clear();
//    d.axes.clear();
//    d.data.resize(isize*jsize, x);
//    d.axes.push_back(isize);
//    d.axes.push_back(jsize);
//    return 0;
//}

//
// find the scalar product
//
template<typename T>
double nrm2(const BOData<T>& data, 
        unsigned int i, 
        unsigned int j, 
        size_t N, 
        double D, 
        double A, 
        double phi, 
        double w, 
        double delT){

    if (N<1)
        return 0.0;

    double scale = 0.0;
    double ssq = 1.0;
    
    size_t isize=data.axes[0];
    size_t jsize=data.axes[1];

    for (unsigned int k = 0; k < N; k++) {
        const double absXk = std::abs(D + A*std::cos(w*i*delT+phi) - data.data[i+isize*(j+jsize*k)]);
        if (absXk != 0.0) {
            if (scale < absXk) {
                ssq = 1.0 + ssq*(scale/absXk)*(scale/absXk);
                scale = absXk;
            } else {
                ssq += (absXk/scale)*(absXk/scale);
            }
        }
    } 
    return scale*std::sqrt(ssq);
}


unsigned int get_signal_index(double signal_frequency, double delT, unsigned int n) {
    double sampling_rate = 1.0/delT;
    double nyquist_frequency = sampling_rate/2.0;
    double del_frequency = sampling_rate/n;

    unsigned int count = 0;
    double frequency = count*del_frequency;

    unsigned int signal_index = 0;

    while (frequency <= nyquist_frequency) {
        std::cerr << count << "   " << frequency;
        if (std::abs(signal_frequency - frequency)/del_frequency < 0.000001) {
            signal_index = count;
            std::cerr << "    *";
        }
        std::cerr << std::endl; 
        ++count;
        frequency = count*del_frequency;
    }

    // need error if signal_index not set.
    return signal_index;
}

template<typename T>
int getAPhiDCFFTW(BOData<T>& data, 
        double freq,    // signal frequency
        unsigned int n,       // samples
        double delT,    // sample interval
        BOData<T>& A, 
        BOData<T>& phi, 
        BOData<T>& DC,
        BOData<T>& corrPhi,
        BOData<T>& snr) {

    size_t isize=data.axes[0];
    size_t jsize=data.axes[1];
    size_t ksize=data.axes[2];

    // Initialize output data
    resetBOData2(A, isize, jsize, static_cast<T>(NAN));
    resetBOData2(phi, isize, jsize, static_cast<T>(NAN));
    resetBOData2(DC, isize, jsize, static_cast<T>(NAN));
    resetBOData2(corrPhi, isize, jsize, static_cast<T>(NAN));
    resetBOData2(snr, isize, jsize, static_cast<T>(NAN));

    unsigned int si = get_signal_index(freq, delT, n); 

    std::cerr << "signal index: " << si << std::endl;

    double* in = reinterpret_cast<double*>(fftw_malloc(sizeof(double)*n)); 
    double* out = reinterpret_cast<double*>(fftw_malloc(sizeof(double)*n)); 
    fftw_plan p = fftw_plan_r2r_1d(n, in, out, FFTW_R2HC, FFTW_MEASURE);

    std::cerr << "Created the FFTW plan" << std::endl;
    double a;
    for (unsigned int i=0; i<data.axes[0]; ++i) {
        for (unsigned int j=0; j<data.axes[1]; ++j) {

            unsigned int idx    = i+isize*j;
            //std::cerr << "idx: " << idx << std::endl;

            for (unsigned int k=0; k<n; ++k) {
                in[k] = data.data[i+isize*(j+jsize*k)];
            }
            fftw_execute(p);
            DC.data[idx] = std::abs(out[0])/n;
            A.data[idx] =  std::sqrt(out[si]*out[si] + out[n-si]*out[n-si])/(n/2);
            phi.data[idx] = std::atan2(out[n-si], out[si]);

            snr.data[idx] = 0.0;
            for (unsigned int k=1; k<(n+1)/2; ++k) {
                a = std::sqrt(out[k]*out[k] + out[n-k]*out[n-k])/(n/2);
                snr.data[idx] += sqrt(a*a); // Non-squared amplitude (modification from fftwfit20110603.h)
            }
            if (n % 2 == 0) { /* N is even */
                a = std::abs(out[n/2])/n;  /* Nyquist freq. */
                snr.data[idx] += sqrt(a*a); // Non-squared amplitude
            }
            snr.data[idx] = sqrt(A.data[idx]*A.data[idx])/(snr.data[idx] - sqrt(A.data[idx]*A.data[idx]));// Non-squared amplitude (modification from fftwfit20110603.h)

            corrPhi.data[idx] = phi.data[idx] + (M_PI/2.0);
            while (corrPhi.data[idx] > (2*M_PI)) {
                corrPhi.data[idx]  -= (2*M_PI);
            }
            while (corrPhi.data[idx] < 0.0) {
                corrPhi.data[idx]  += (2*M_PI);
            }
            //
            // Do error estimates
            //
        }
    }
    fftw_destroy_plan(p);
    if (in) 
        fftw_free(in);
    if (out)
        fftw_free(out);

    return 0;
}
#endif 
