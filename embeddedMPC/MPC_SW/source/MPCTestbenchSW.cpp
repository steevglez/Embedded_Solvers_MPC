#include <iostream>
#include <fstream>
#include "specs.h"
#include <math.h>
#include "mpcSW.h"
#include "tictoc.h"
#include <iomanip>


using namespace std;
#define DISPLAY

int main(int argc, char *argv[])
{
    cout << std::setprecision(15);
    cout << "MPC testbench" << endl;

	elem x[N_SYS];
	elem yref[P_SYS];
	elem xmin[N_SYS];
	elem xmax[N_SYS];
	elem umin[M_SYS];
	elem umax[M_SYS];
	elem Acal[N_SYS*N_QP][N_SYS];
	elem AcalOmgOcal[N_SYS][N_QP];
	elem H[N_QP *N_QP];
	elem Mx[I_QP *N_QP];
	elem L_invLast[(N_SYS+P_SYS)];
	elem expectedResult[P_SYS];
	elem u[P_SYS];
	int iterPDIP = 20;
	int iterLS = 30;
	int errors = 0;

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

    int nSamples = 0;
    int samplesI = 0;
    int samplesN = 0;


    // tic toc timer
    tictoc MPCTimer(nSamples);
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

    // load xmin
    for (int r=0; r<N_SYS; r++){
    	samples >> xmin[r];
    }
    // load xmax
    for (int r=0; r<N_SYS; r++){
    	samples >> xmax[r];
    }
    // load umin
    for (int r=0; r<M_SYS; r++){
    	samples >> umin[r];
    }
    // load umax
    for (int r=0; r<M_SYS; r++){
    	samples >> umax[r];
    }
    // load Acal
    for (int r=0; r<N_SYS*N_QP; r++){
    	for (int c=0; c<N_SYS; c++){
    		samples >> Acal[r][c];
    	}
    }
    // load AcalOmgOcal
    for (int r=0; r<N_SYS; r++){
    	for (int c=0; c<N_QP; c++){
    		samples >> AcalOmgOcal[r][c];
    	}
    }
    // load H
    for (int r=0; r<N_QP; r++){
    	for (int c=0; c<N_QP; c++){
    		samples >> H[r*N_QP +c];
    	}
    }
    // load Mx
	for (int r=0; r<I_QP; r++){
		for (int c=0; c<N_QP; c++){
			samples >> Mx[r*N_QP +c];
		}
    }
    // load L_last
	for (int r=0; r<N_SYS+P_SYS; r++){
		samples >> L_invLast[r];
    }

    for (int sample=0; sample<nSamples; sample++){
    	// exit if end of file
        if (samples.peek() == EOF) break;
        // load x
        for (int i=0; i<N_SYS; i++){
            samples >> x[i];
        }
        // load yref
        for (int i=0; i<P_SYS; i++){
            samples >> yref[i];
        }
        // load expectedResult
        for (int i=0; i<P_SYS; i++){
            samples >> expectedResult[i];
        }

        MPCTimer.tic();
        mpcSW(x, yref, xmin, xmax, umin, umax, Acal, AcalOmgOcal, H, Mx, L_invLast, iterPDIP, iterLS, tol, u);
        MPCTimer.toc();

        //cout << "sample number : " << sample << endl;
        for (int i=0; i<P_SYS; i++){
        	elem error = fabs(expectedResult[i] -u[i]);
            if ((error > threshold) | (error != error)){
                errors++;
#ifdef DISPLAY
                cout << "sample number : " << sample << " " << i << endl;
                cout << expectedResult[i] << "\t";
                cout << u[i] << "\t";
                cout << error;
                cout << "\n";
#endif

            }
        }

    }

    samples.close();

    // print important data
    cout << "Finished processing " << nSamples << " samples" << endl;
    cout << "Number of differences between expected and calculated:\t" << errors << endl;
    cout << "Threshold: " << threshold << endl;
    std::cout << std::endl;

    MPCTimer.rmWarmup(1000);
    double sw_cycles = MPCTimer.mean();

    std::cout << "Mean time running application in software: \t\t" << sw_cycles / 1e3 << "\xC2\xB5s" << std::endl;
    std::cout << "Standard deviation running application in software: \t" << MPCTimer.stdev() / 1e3 << "\xC2\xB5s" << std::endl;
    std::cout << "Max running application in software: \t\t\t" << MPCTimer.max() / 1e3 << "\xC2\xB5s" << std::endl;
    std::cout << "Min running application in software: \t\t\t" << MPCTimer.min() / 1e3 << "\xC2\xB5s" << std::endl;

    return errors;
}
