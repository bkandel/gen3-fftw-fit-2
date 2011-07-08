#include <valarray>
#include <string>
#include <gen3pp/bodata.h>
#include <gen3pp/boutil.h>
#include <gen3pp/ccfitsfile.h>
#include <gen3pp/fftwfit.h>
#include "paramlist.h"

using namespace std;

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
    

    BOData<double> sigma;
    double preset_sigma;
    unsigned int err = pl.getValue("RF_PRESET_SIGMA", preset_sigma);
    if (err == 0) 
    {
        cerr << "Processing Preset Sigma ...";
        sigma.data = darkmean.data;
        sigma.axes = darkmean.axes;
        string filename = "sigmasigma3.txt";
        write2file2(filename, sigma);
        cerr << " done." << endl;
    } 
    else 
    {
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
    
    
	double freq;
    double del_t;
    double nOfData;
    const char *wisdom_filename = "wisdom.txt"; 
    
    pl.getValue("RF_DEL_TIME", del_t);
    pl.getValue("RF_MODULATED_FREQ", freq);
    pl.getValue("RF_NUM_OF_DATA_POINTS", nOfData);
    
    cerr << "Creating FFTW Plan..." << endl;
    
    plan_fftw(nOfData, wisdom_filename); 
    cerr << "FFTW Plan creation done." << endl; 
    
    return 0;
}
