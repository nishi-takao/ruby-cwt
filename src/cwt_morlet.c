//
// Continuous wavelet transform using Morlet wavelet for NArray
//
// Time-stamp: <2013-02-12 18:28:18 zophos>
//
// Copyright (c) 2013 NISHI, Takao <zophos@ni.aist.go.jp> AIST
// All rights reserved.
//
// This is free software with ABSOLUTELY NO WARRANTY.
//
// You can redistribute it and/or modify it
// under the terms of GPLv2 or later.
//
//
// originate from "Rwave: Time-Frequency analysis of 1-D signals"
// http://cran.r-project.org/web/packages/Rwave/
//
/****************************************************************
*               (c) Copyright  1997                             *
*                          by                                   *
*      Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   *
*                  Princeton University                         *
*                  All right reserved                           *
****************************************************************/

#include <ruby.h>
#include <narray.h>

#include <math.h>
#include <stdlib.h>

static const double TWOPI = 6.28318530717959;


static int find2power(int n)
{
    long m, m2;

    m = 0;
    m2 = 1<<m; /* 2 to the power of m */
    while (m2-n < 0) {
        m++;
        m2 <<= 1; /* m2 = m2*2 */
    }
    return(m);
}

/******************************************
*  FFT function (from Numerical Recipes)  *
*  (the length of the FFTed signal HAS    *
*  be a power of 2)                       *
*******************************************/
static void four1(double data[], int nn, int isign)
{
    int n,mmax,m,j,istep,i;
    double wtemp,wr,wpr,wpi,wi,theta;
    double tempr,tempi;
    
    n=nn << 1;
    j=1;
    for (i=1;i<n;i+=2) {
        if (j > i) {
            // SWAP(data[j],data[i]);
            tempr=data[j];
            data[j]=data[i];
            data[i]=tempr;
        
            //SWAP(data[j+1],data[i+1]);
            tempr=data[j+1];
            data[j+1]=data[i+1];
            data[i+1]=tempr;
        }

        m=n >> 1;
        while (m >= 2 && j > m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }

    mmax=2;
    while (n > mmax) {
        istep=2*mmax;
        theta=TWOPI/(isign*mmax);
        wtemp=sin(0.5*theta);
        wpr = -2.0*wtemp*wtemp;
        wpi=sin(theta);
        wr=1.0;
        wi=0.0;
        for (m=1;m<mmax;m+=2) {
            for (i=m;i<=n;i+=istep) {
                j=i+mmax;
                tempr=wr*data[j]-wi*data[j+1];
                tempi=wr*data[j+1]+wi*data[j];
                data[j]=data[i]-tempr;
                data[j+1]=data[i+1]-tempi;
                data[i] += tempr;
                data[i+1] += tempi;
            }
            wr=(wtemp=wr)*wpr-wi*wpi+wr;
            wi=wi*wpr+wtemp*wpi+wi;
        }
        mmax=istep;
    }
}


/* Fast Fourier transform (from numerical recipes routine)
   -------------------------------------------------------
*/
static void double_fft(double *Or,double *Oi,double *Ir,double *Ii,
                       double *tmp,
                       int isize,int newsize,int isign)
{
    int i;
    double mag=1.0;

    for(i = 0; i < isize; i++) {
        tmp[2 * i] = Ir[i];
        tmp[2 * i + 1] = Ii[i];
    }
    four1(tmp-1,newsize,isign);  
    
    
    if(isign == -1)
        mag=1.0/newsize;
    else
        mag=1.0;
    
    for(i = 0; i < isize; i++) {
        Or[i] = tmp[2 * i] * mag;
        Oi[i] = tmp[2 * i + 1] * mag;
    }
}

/***************************************************************
*  Function: multi:
*  ---------
*     Multiplication of 2 vectors in the Fourier domain; the
*     first one is complex-valued, as well as the output; the
*     second one is real-valued.
*
*   Ri1, Ii1: first input vector.
*   Ri2: second input vector (real).
*   Or, Oi: output vector.
*   isize: length of the vectors.
***************************************************************/

static void multi(double *Ri1, double *Ii1, double *Ri2, double *Or,
                  double *Oi, int isize)
{
    int i;
    
    for(i = 0; i < isize; i++) {
        Or[i] = Ri1[i] * Ri2[i];
        Oi[i] = Ii1[i] * Ri2[i];
    }
    return;
}


/***************************************************************
*  Function: morlet_frequency:
*  ---------
*     Generates a Morlet wavelet in the frequency domain.
*     The wavelet is centered at the origin, and normalized so
*     that psi(0) =1/sqrt(2 pi)
*
*   w: wavelet
*   scale: scale of the wavelet
*   isize: window size
*   cf: central frequency of the wavelet
***************************************************************/

static void morlet_frequency(double cf,double scale,double *w,int isize)
{
    double tmp;
    int i;
    
    for(i = 0; i < isize; i++) {
        tmp = (double)(scale * i * TWOPI/isize - cf);
        tmp = -(tmp * tmp)/2;
        w[i] = exp(tmp);
    }
    return;
}


/*****************************************************************
*  function:  Scwt_morlet
*    Continuous wavelet transform :
*
*   input: (a priori complex-valued) input signal
*   Ri1, Ii1: Fourier transform of input signal (real and
*      imaginary parts).
*   Ri2: Real part of Fourier transform of Morlet wavelet
*   Oreal,Oimage: real and imaginary parts of CWT
*   inputsize: signal size
*   nboctave: number of scales (powers of 2)
*   nvoice: number of scales between 2 consecutive powers of 2
*   centerfrequency: centralfrequency of Morlet wavelet
******************************************************************/

void Scwt_morlet(double *Rinput,double *Iinput,
                 double *Oreal,double *Oimage,
                 int nboctave,int nbvoice,int inputsize,
                 double centerfrequency)
{
    int i, j,fftsize;
    double a;
    double *Ri2, *Ri1, *Ii1, *Ii, *Ri, *FFTBuf, *WorkBuffers;
    

    fftsize=(1 << find2power(inputsize));

    //
    // allocate all buffers to avoid memory leak.
    //
    WorkBuffers=(double *)ruby_xcalloc(inputsize*5+fftsize*2,sizeof(double));

    Ri2=WorkBuffers;
    Ri1=Ri2+inputsize;
    Ii1=Ri1+inputsize;
    Ri=Ii1+inputsize;
    Ii=Ri+inputsize;
    FFTBuf=Ii+inputsize;

    for(i = 0; i < inputsize; i++) {
        Ri[i] = (double)Rinput[i]; 
        Ii[i] = (double)Iinput[i];
    }

    double_fft(Ri1,Ii1,Ri,Ii,FFTBuf,inputsize,fftsize,-1);
    

    for(i = 1; i <= nboctave; i++) {
        for(j=0; j < nbvoice; j++) {
            a = (double)(pow((double)2,(double)(i+j/((double)nbvoice))));
            morlet_frequency(centerfrequency,a,Ri2,inputsize); 
            multi(Ri1,Ii1,Ri2,Oreal,Oimage,inputsize);
            double_fft(Oreal,Oimage,Oreal,Oimage,FFTBuf,inputsize,fftsize,1); 
            Oreal = Oreal + inputsize;
            Oimage = Oimage + inputsize;  
        }
    }

    ruby_xfree(WorkBuffers);
}


VALUE rb_na_cwt_morlet(VALUE self,
                       VALUE rbo_noctave,
                       VALUE rbo_nvoice,
                       VALUE rbo_w0)
{

    
    struct NARRAY *na_input;
    GetNArray(self,na_input);

    //
    // check self.total
    //
    int size=na_input->total;

    double oct=log2((double)size);
    if(oct!=(double)((int)oct))
        rb_raise(rb_eRuntimeError,"total size must be equal to 2**n.");

    //
    // check noctave
    //
    int noctave=FIX2INT(rbo_noctave);
    oct=log2((double)noctave);
    if(oct!=(double)((int)oct))
        rb_raise(rb_eRuntimeError,"noctave must be equal to 2**n.");


    //
    // decomposite complex array to real and image array
    //
    VALUE rbo_input;
    if(na_input->type==NA_DCOMPLEX)
        rbo_input=self;
    else
        rbo_input=na_cast_object(self,NA_DCOMPLEX);

    VALUE rbo_input_real=rb_funcall(rbo_input,
                                    rb_intern("real"),
                                    0);
    VALUE rbo_input_imag=rb_funcall(rbo_input,
                                    rb_intern("imag"),
                                    0);


    //
    // prepare output buffer
    //
    int nvoice=FIX2INT(rbo_nvoice);
    VALUE pp=INT2FIX(noctave*nvoice);

    VALUE cNMatrix=rb_const_get(rb_cObject,rb_intern("NMatrix"));
    VALUE rbo_output_real=rb_funcall(cNMatrix,
                                     rb_intern("dfloat"),
                                     2,
                                     INT2FIX(size),
                                     pp);
    VALUE rbo_output_imag=rb_funcall(cNMatrix,
                                     rb_intern("dfloat"),
                                     2,
                                     INT2FIX(size),
                                     pp);
    
    struct NARRAY *na_input_real;
    GetNArray(rbo_input_real,na_input_real);
    struct NARRAY *na_input_imag;
    GetNArray(rbo_input_imag,na_input_imag);

    struct NARRAY *na_output_real;
    GetNArray(rbo_output_real,na_output_real);
    struct NARRAY *na_output_imag;
    GetNArray(rbo_output_imag,na_output_imag);
    
    
    Scwt_morlet((double *)na_input_real->ptr,
                (double *)na_input_imag->ptr,
                (double *)na_output_real->ptr,
                (double *)na_output_imag->ptr,
                noctave,
                nvoice,
                size,
                NUM2DBL(rbo_w0));

    //
    // composite complex array from real and image array
    //
    VALUE rbo_output=na_cast_object(rbo_output_real,NA_DCOMPLEX);
    
    rb_funcall(rbo_output,
               rb_intern("imag="),
               1,
               rbo_output_imag);

    return rbo_output;
}

void Init_cwt_morlet()
{
    rb_require("narray");
    VALUE klass=rb_const_get(rb_cObject, rb_intern("NVector"));
    rb_define_method(klass,"cwt_morlet",rb_na_cwt_morlet,3);
}
