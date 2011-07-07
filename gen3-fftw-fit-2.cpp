#include <valarray>
#include <string>
#include <gen3pp/bodata.h>
#include <gen3pp/boutil.h>
#include <gen3pp/ccfitsfile.h>
#include <gen3pp/fftwfit.h>
#include "paramlist.h"

using namespace std;

const int sigFig = 10; 

int main(int argc, char* argv[]) {

    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <filename.cfg>" << endl;
        return 0;
    }

    ParamList pl; 
    cerr << "Reading configuration...";
    pl.parse(argv[1]);
    cerr << " done." << endl;

    BOData<double> darkmean;

    {
        cerr << "Processing Dark frame...";
        string filename;
        pl.getValue("RF_DARK_FILE", filename);
        CCFITSfile fitsfile;
        BOData<double> dat;
        fitsfile.read(filename, dat);
        print(dat);

        getMean(dat, 3, darkmean); 
        filename = "darkmean3.txt";
        write2file2(filename, darkmean);
        cerr << " done." << endl;
    }

    BOData<double> sigma;
    double preset_sigma;
    unsigned int err = pl.getValue("RF_PRESET_SIGMA", preset_sigma);
    if (err == 0) {
        cerr << "Processing Preset Sigma ...";
        sigma.data = darkmean.data;
        sigma.axes = darkmean.axes;
        string filename = "sigmasigma3.txt";
        write2file2(filename, sigma);
        cerr << " done." << endl;
    } else {
        cerr << "Processing Sigma frame...";
        string filename;
        pl.getValue("RF_SIGMA_FILE", filename);
        CCFITSfile fitsfile;
        BOData<double> sigmaData;
        fitsfile.read(filename, sigmaData);
        print(sigmaData);
        getSD(sigmaData, 3, sigma);
        filename = "sigmasigma3.txt";
        write2file2(filename, sigma);
        cerr << " done." << endl;
    }
        
    BOData<double> dat;
    {
        cerr << "Processing Data frame...";
        string filename;
        pl.getValue("RF_DATA_FILE", filename);
        CCFITSfile fitsfile;
        fitsfile.read(filename, dat);
        print(dat);
        cerr << " done." << endl;

        cerr << "Subtracting Dark frame...";
        subtract(darkmean, dat);
        cerr << " done." << endl;
    }

    BOData<double> A;
    BOData<double> phi;
    BOData<double> DC;
    BOData<double> corrPhi;
    BOData<double> snr;

    {
        double freq;
        double del_t;
        double nOfData;
        double absErr;
        double relErr;
        unsigned int maxIter;
        pl.getValue("RF_DEL_TIME", del_t);
        pl.getValue("RF_MODULATED_FREQ", freq);
        pl.getValue("RF_NUM_OF_DATA_POINTS", nOfData);

        cerr << "Fitting for A, Phi and DC..." << endl;
        //getAPhiDC(dat, sigma, del_t, freq, nOfData, absErr, relErr, maxIter, A, phi, DC, iters, Status, delDC, delA, delPhi, chi2ByDOF, initA, initPhi, initDC);
        getAPhiDCFFTW(dat, freq, nOfData, del_t, A, phi, DC, corrPhi, snr);      // signal frequency
        cerr << " done." << endl;
    }

    string filename = "corrPhi.gim";
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