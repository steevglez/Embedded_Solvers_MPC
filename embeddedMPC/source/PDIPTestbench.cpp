#include <iostream>
#include <fstream>
#include "specs.h"
#include "pdipHW.h"
#include <iomanip>


using namespace std;
#define DISPLAY

int main(int argc, char *argv[])
{
    cout << std::setprecision(15);
    cout << "PDIP testbench" << endl;
    elem H[N_QP][N_QP];
    elem h[N_QP];
    elem Mx[I_QP][N_QP];
    elem cx[I_QP];
    elem expectedResult[N_QP];
    elem x[N_QP];
    elem zko[N_QP];

    int errors = 0;
    for (int i=0; i<N_QP; i++){
        zko[i] = 0.0;
    }

    if (argc!=2){
        cerr << "Must specify .txt\n";
        return 1;
    }

    ifstream samples(argv[1], ios::in);
    //check to see that the file was opened correctly:
    if (!samples.is_open()) {
        cerr << "There was a problem opening the input file: ";
        cerr << argv[1] << endl;
        return 1;//exit or do additional error checking
    }

    // threshold for the difference between the calculated results and the correct ones.
    elem threshold = 1e-4;
    // tolerance for the linear solver in PDIP algorithm
    elem tol = 1e-9;
    // PDIP iterations
    int iterPDIP = 20;
    // MINRES iterations
    int iterMINRES = 30;

    int nSamples = 0;
    int samplesI = 0;
    int samplesN = 0;

    // number of samples to read from samples file
    // large numbers can make the test bench really slow
    int nSamplestb = 5000;


    string line;
    // get number N, I and number of samples
    getline(samples, line);
    sscanf(line.c_str(), "N:\t%d\tI: %d\tnSamples:\t%d", &samplesN, &samplesI, &nSamples);

    // check if N and I match with N_QP and I_QP from specs.h
    if (N_QP != samplesN){
        cerr << "N_QP from " << argv[1] << " (" << samplesN << ") does not match N_QP from specs.h (";
        cerr << N_QP << ")" << endl;
        return 1;
    }
    if (I_QP != samplesI){
        cerr << "I_QP from " << argv[1] << " (" << samplesI << ")does not match I_QP from specs.h (";
        cerr << I_QP << ")" << endl;
        return 1;
    }


    // load H
    for (int r=0; r<N_QP; r++){
    	for (int c=0; c<N_QP; c++){
    		samples >> H[r][c];
    	}
    }
    // load Mx
	for (int r=0; r<I_QP; r++){
		for (int c=0; c<N_QP; c++){
			samples >> Mx[r][c];
		}
    }


    for (int sample=0; sample<nSamplestb; sample++){
    	// exit if end of file
        if (samples.peek() == EOF) break;
        // load h
        for (int i=0; i<N_QP; i++){
            samples >> h[i];
        }
        // load cx
        for (int i=0; i<I_QP; i++){
            samples >> cx[i];
        }
        // load expectedResult
        for (int i=0; i<N_QP; i++){
            samples >> expectedResult[i];
        }

        pdipHW(H, h, Mx, cx, zko, x, iterPDIP, iterMINRES, tol);

        // reset zko (initial linear system solution)
        for (int i=0; i<N_QP; i++){
            zko[i] = 0.0;
        }

        //cout << "sample number : " << sample << endl;
        for (int i=0; i<N_QP; i++){
        	elem error = fabs(expectedResult[i] -x[i]);
        	//cout << error<< "\n";
            if ((error > threshold) | (error != error)){
                errors++;
#ifdef DISPLAY
                cout << "sample number : " << sample << " " << i << endl;
                cout << expectedResult[i] << "\t";
                cout << x[i] << "\t";
                cout << error;
                cout << "\n";
#endif

            }
        }
        
    }

    samples.close();

    // print important data
    cout << "Finished " << nSamplestb << " QP problems" << endl;
    cout << "Number of differences between expected and calculated " << errors << endl;
    cout << "Threshold: " << threshold << endl;
    std::cout << std::endl;

    return errors;
}
