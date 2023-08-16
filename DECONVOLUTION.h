#ifndef DECONVOLUTION_H
#define DECONVOLUTION_H
#include <fftw3.h>
#include "constants.h"
#include <TMath.h>
#define Re 0
#define Im 1
class DECONVOLUTION
{
    public:

    void SetZl(Int_t zero_level)
    {
        zl = zero_level;
    }
    void SetSignal(Short_t *in)
    {
        for (int i = 0; i < SNAPSHOT_LENGTH; i++) 
        {
            inSignal[i][Re] = (double)in[i]-zl;
            inSignal[i][Im] = 0;
        }
    }

    template<typename T>
    void GetSignal(T *in)
    {

        for (int i = 0; i < SNAPSHOT_LENGTH; i++) 
        {
            in[i] = 0;
            for (int j = 0; j <= i; j++) in[i] += (Short_t)modcomplex(inSignal[j]);
        }    
        Float_t zl_local = 0;
        for (int i = 0; i < ZERO_LEVEL_POS; i++) zl_local += in[i]/ZERO_LEVEL_POS;
        for (int i = 0; i < SNAPSHOT_LENGTH; i++) in[i]=in[i]-zl_local+zl;

    }

    void SetResponse()
    {
        float tau = 235;
        for (int i =0; i < SNAPSHOT_LENGTH; i++)
        {
            inResponse[i][Re]= exp(-((float)i)/tau);
            inResponse[i][Im]= 0;
        }
        fft(inResponse,transResponse);
    }


    // void PrintResponseReIm()
    // {
    //     SetResponse();
    //     cout << "\n";
    //     fft(inResponse,transResponse);
    //     ifft(transResponse,inResponse);
    //     for (int i = 0; i < SNAPSHOT_LENGTH; i++) 
    //         cout << ((inResponse[i][Re] > 0) - (inResponse[i][Re] < 0) )*
    //         sqrt(pow(inResponse[i][Re],2) + pow(inResponse[i][Im],2)) << "\t";       
    // }

    void deconvolution()
    {
        fft(inSignal,transSignal);
        complex_divide(transSignal,transResponse);
        ifft(transSignal,inSignal);
    }   

    void Reset()
    {
        for (int i = 0; i < SNAPSHOT_LENGTH; i++)
        {
            inSignal[i][Re] = 0;
            transSignal[i][Re] = 0;
            // inResponse[i][Re] = 0;
            // transResponse[i][Re] = 0;   
            inSignal[i][Im] = 0;
            transSignal[i][Im] = 0;
            // inResponse[i][Im] = 0;
            // transResponse[i][Im] = 0;           
        }
        zl = 0;
    }
    private:
    fftw_complex inSignal[SNAPSHOT_LENGTH];
    fftw_complex transSignal[SNAPSHOT_LENGTH];
    fftw_complex inResponse[SNAPSHOT_LENGTH];
    fftw_complex transResponse[SNAPSHOT_LENGTH];
    Short_t zl=0;

    void fft(fftw_complex *in, fftw_complex *out)
    {
        fftw_plan plan = fftw_plan_dft_1d(SNAPSHOT_LENGTH,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);
        fftw_cleanup();
        for (int i = 0; i < SNAPSHOT_LENGTH; i++)
        {
            // out[i][Re] /=sqrt(SNAPSHOT_LENGTH);
            // out[i][Im] /=sqrt(SNAPSHOT_LENGTH);
        }
    }

    void ifft(fftw_complex *in, fftw_complex *out)
    {
        fftw_plan plan = fftw_plan_dft_1d(SNAPSHOT_LENGTH,in,out,FFTW_BACKWARD,FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);
        fftw_cleanup();
        for (int i = 0; i < SNAPSHOT_LENGTH; i++)
        {
            out[i][Re] /=SNAPSHOT_LENGTH;
            out[i][Im] /=SNAPSHOT_LENGTH;
        }
    }
    void complex_divide(fftw_complex *in, fftw_complex *out)
    {
        fftw_complex temp;
        for (int i = 0; i < SNAPSHOT_LENGTH; i++)
        {

            temp[Re] = (in[i][Re]*out[i][Re]+in[i][Im]*out[i][Im])
            /(out[i][Re]*out[i][Re]+out[i][Im]*out[i][Im]);
            temp[Im] = (in[i][Im]*out[i][Re]-in[i][Re]*out[i][Im])
            /(out[i][Re]*out[i][Re]+out[i][Im]*out[i][Im]);


            in[i][Re] = temp[Re];
            in[i][Im] = temp[Im];            
        }
    }

    double modcomplex(fftw_complex a)
    {
        return ((a[Re] > 0) - (a[Re] < 0) )*
            sqrt(pow(a[Re],2) + pow(a[Im],2));
    }

};

#endif DECONVOLUTION_H