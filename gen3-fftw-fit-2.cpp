#include <valarray>
#include <string>
#include <gen3pp/bodata.h>
#include <gen3pp/boutil.h>
#include <gen3pp/ccfitsfile.h>
#include <gen3pp/fftw2fit.h>
#include "paramlist.h"
#include <time.h>

using namespace std;

const int sigFig = 10; 
std::string darkfile = "darkmean3.gim"; 

int main(int argc, char* argv[]) {

    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <filename.cfg>" << endl;
        return 0;
    }
    
    

    ParamList pl; 
    cerr << "Reading configuration...";
    pl.parse(argv[1]);
    cerr << " done." << endl;


    BOData<double> dat;
    cerr << "Processing Data frame...";
    string filename;
    pl.getValue("DATA_FILE", filename);
    CCFITSfile fitsfile;
    fitsfile.read(filename, dat);
    print(dat);
    cerr << " done." << endl;

	BOData<double> darkmean; 
	cerr << "Reading Dark frame..." << endl; 
	readgim(darkfile, darkmean);
    cerr << "Subtracting Dark frame...";
    subtract(darkmean, dat);
    cerr << " done." << endl;
    
    
    BOData<double> A;
    BOData<double> phi;
    BOData<double> DC;
    BOData<double> corrPhi;
    BOData<double> snr;
    const char *wisdom_filename = "wisdom.txt"; 
    double freq;
    double del_t;
    double nOfData;
    double absErr;
    double relErr;
    unsigned int maxIter;
    pl.getValue("SAMPLING_INTERVAL", del_t);
    pl.getValue("XCORR_FREQUENCY", freq);
    pl.getValue("NUM_DATA_POINTS", nOfData);

    cerr << "Fitting for A, Phi and DC..." << endl;
        //getAPhiDC(dat, sigma, del_t, freq, nOfData, absErr, relErr, maxIter, A, phi, DC, iters, Status, delDC, delA, delPhi, chi2ByDOF, initA, initPhi, initDC);
    getAPhiDCFFTW2(dat, freq, nOfData, del_t, A, phi, DC, corrPhi, snr, wisdom_filename);      
    cerr << " done." << endl;



    filename = "corrPhi.gim";
    write2gim(filename.c_str(), corrPhi, sigFig);

    filename = "A.gim";
    write2gim(filename.c_str(), A, sigFig);

    filename = "phi.gim";
    write2gim(filename.c_str(), phi, sigFig);

    filename = "DC.gim";
    write2gim(filename.c_str(), DC, sigFig);

    filename = "snr.gim";
    write2gim(filename.c_str(), snr, sigFig);

    return 0;
}
