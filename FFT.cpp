// FFT.cpp
// Program for computing the Fast Fourier Transform (FFT) of a 1D signal
// Edwin Cordon

#include <iostream>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>

#include "Complex.h"

#define PI	    3.14159265359
#define FFT     1

#define debug   false

typedef unsigned int uint;

using namespace std;

void DFT(Complex X[], double x[], uint N);
void FFT1D(Complex X[], double x[], uint N);

uint getPaddedLength(uint length); 
void createPaddedx(double padx[], double x[], uint paddedLength);

int trace(bool on, const char *format, ...);
int readInput(const char *fname, double x[], uint N);
int saveFFT(const char *sname, Complex X[], uint N);

int main(int argc, char *argv[]){

    char *sname;
    int i = 1;
    uint samples;
    double *x;
    
    if (argc != 6){
        fprintf(stderr, "./FFT -i <samples> <input file> -o <output file>\n");
        return -1;
    }

    while (i < argc){

        if (argc - i >= 3){
        
            if ((0 == strcmp(argv[i], "-i"))){

                int status;

                trace(debug, "argv[%d] = %s, argv[%d] = %s, argv[%d] = %s\n",
                        i, argv[i], i+1, argv[i+1], i+2, argv[i+2]);

                sscanf(argv[i+1], "%u", &samples); 
                x = new double[samples];
                
                status = readInput(argv[i+2], x, samples);

                if (status != 0)
                    return status;

                i += 3;
            }   
        }

        if (argc - i >= 2){

            if ((0 == strcmp(argv[i], "-o"))){

                trace(debug, "argv[%d] = %s, argv[%d] = %s\n",
                        i, argv[i], i+1, argv[i+1]);
                
                // This can be done better with C++, but
                // right now I'm just going with what I know.

                sname = (char *)malloc(strlen(argv[i+1]));
                strcpy(sname, argv[i+1]);
                
                trace(debug, "sname = %s\n", sname);

                i += 2;
            }
        }

    }

#if FFT // Calculate Fast Fourier Transform if FFT flag is set

    uint paddedLength = getPaddedLength(samples);

    double *padx = new double [paddedLength]; 
	Complex *X = new Complex[paddedLength];

    createPaddedx(padx, x, paddedLength);

    FFT1D(X, padx, paddedLength);

    saveFFT(sname, X, paddedLength);

    delete [] padx;
    delete [] X;

#else // Otherwise calculate the Discrete Fourier Transform

    Complex X[samples];

	DFT(X, x, samples);

	saveFFT(sname, X, samples);
#endif

    delete [] x;
    free(sname);

	return 0;
}

void DFT(Complex X[], double x[], uint N){

	for(uint k = 0; k < N; k++){
		for(uint n = 0; n < N; n++){

			Complex W(1.0, -2*PI*n*k/N);
			X[k] += W*x[n];
		}
	}
}

void FFT1D(Complex X[], double x[], uint N){

    trace(debug, "\nStarting %u-point FFT\n", N);

    if(N == 1){
		
        X[0] = x[0];
	}
    else if (N%2 == 0){
               
        Complex *EVEN = new Complex[N/2];
		Complex *ODD = new Complex[N/2];

    	double *even = new double[N/2];
	    double *odd = new double[N/2];

		// spilt even and odd samples from input signal
        for(uint n = 0; n < N/2; n++){

	    	even[n] = x[2*n];
		    odd[n] = x[2*n+1];
    
            trace(debug, "x[%u] = %lf, x[%u] = %lf\n", 
                    2*n, x[2*n], 2*n+1, x[2*n+1]);

            trace(debug, "even[%u] = %lf, odd[%u] = %lf\n",
                    n, even[n], n, odd[n]);
        }

        // get even and odd N/2-point FFTs
		FFT1D(EVEN, even, N/2);
    	FFT1D(ODD, odd, N/2);

	    delete [] even;
		delete [] odd;

        trace(debug, "\nResult for %u-point FFT\n", N);

        // get N-point FFT
		for(uint k = 0; k < N/2; k++){
    
            Complex Wl(1.0, -2*PI*k/N);
            Complex Wh(1.0, -2*PI*(k)/N); 

            X[k] = EVEN[k] + Wl*ODD[k];
            X[k+N/2] = EVEN[k] - Wh*ODD[k];

            trace(debug, "|X[%u]| = |(%lf + j(%lf)) + (%lf + j(%lf))(%lf - j(%lf))| = %lf\n", 
                    k,
                    EVEN[k].getReal(), EVEN[k].getImaginary(), 
                    Wl.getReal(), Wl.getImaginary(), 
                    ODD[k].getReal(), ODD[k].getImaginary(), 
                    X[k].getMagnitude());

            trace(debug, "|X[%u]| = |(%lf + j(%lf)) + (%lf - j(%lf))(%lf - j(%lf))| = %lf\n", 
                    k+N/2,
                    EVEN[k].getReal(), EVEN[k].getImaginary(),
                    Wh.getReal(), Wh.getImaginary(), 
                    ODD[k].getReal(), ODD[k].getImaginary(),
                    X[k+N/2].getMagnitude());
        }
    
        delete [] EVEN;
		delete [] ODD;    
    }
}

uint getPaddedLength(uint length){

    uint i, mask, numBits;
    uint msb = 0, asserts = 0;

    // Look for most significant bit of the input length
    // to determine the proper padding.
    
    numBits = sizeof(uint)*8;
    mask = (1 << (numBits - 1));

    for (i = 0; i < numBits; i++){
        
        if (((length << i) & mask)){
            
            if (msb == 0)
                msb = numBits - i - 1;

            asserts++;
        }
    }

    if (asserts <= 1)
        return length;
    else
        return (1 << (msb + 1));
}

void createPaddedx(double padx[], double x[], uint paddedLength){

    for(uint n = 0; n < paddedLength; n++){

        if (n < paddedLength)
            padx[n] = x[n];
        else
            padx[n] = 0.0;
    }
}

int trace(bool on, const char *format, ...){

    if (on){

        va_list arg;
        int status;

        va_start(arg, format);
        status = vfprintf(stdout, format, arg);
        va_end(arg);

        return status;
    }
    else
        return 0;
}

int readInput(const char *fname, double x[], uint N){

    FILE *fp;

    try{

        if (!(fp = fopen(fname, "r")))
            throw -1;
    }
    catch (int e){
        
        if (e == -1)
            return e;
    }

    for (uint n = 0; n < N; n++)
        fscanf(fp, "%lf", &x[n]);
    
    fclose(fp);

    return 0;
}

int saveFFT(const char *sname, Complex X[], uint N){

	FILE *fp;

	try{

		if (!(fp = fopen(sname, "w")))
			throw -1;
	}
	catch (int e){
		
		if (e == -1)
			return e;
	}

	for (uint k = 0; k < N; k++)
		fprintf(fp, "%lf ", X[k].getMagnitude());

    fprintf(fp, "\n");

    fclose(fp);

    return 0;
}

