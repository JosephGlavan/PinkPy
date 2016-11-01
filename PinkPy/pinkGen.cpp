// Joseph Glavan j.glavan4@gmail.com
// Running xxx's code through the grinder to make it do what I want


// Not sure which of these are necessary
//#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "pinkGen.h"
//#include <allegro.h>
//#include <winalleg.h>
//#include "gamma.h"
//#include "debug.h"
//#include "ini.h"

//#define M_PI	3.14159

// Changing dimension declarations from #define's to simple int (so that I can parameterize them)
//#define NX			256

//#define N  NX // I don't think this gets used
//#define NY NX
//#define NZ 64
#define NZ_STEP 1
#define NZ_MULT 1
//#define NSIZE (NX * NY * NZ)


/*
int NX;// = 256;
int NY;// = 256;
int NZ;// = 16;
int NSIZE;// = (NX * NY * NZ);
int SeedVal; // = 0;
*/


// These don't appear necessary for me
/*
#define TSLICE_64	33.3333
#define TSLICE_512	16.6666
#define TSLICE		TSLICE_64

int Gamma = 0;
int MaxX = 1024;
int MaxY = 768;
int ColorDepth = 32;
int save = FALSE;
int MaxScale = 63;              // We only have 64 shade when using a pallete    
int mode  = NX;
int gray  = 1;
float gain = 1.0;
float dc   = 0.5;
float mdc   = 0.0;
int WHITE,BLACK,GRAY1,GRAY2,GRAY3,GRAY4,GRAY5,GRAY6,GRAY7,GRAY8,GRAY9;
int showHelp = TRUE;
int showHeader = FALSE;
int build = FALSE;
float popStrength = 1.0;
float popDensity = 0;
int ErrorFlag = 0;
*/

//long TM0,TM1,TM2;
//#define PROFILE(arg1) { TM2 = GetTickCount(); Debug("Time: %6.3f : %s took %6.3f\n",(TM1-TM0)/1000.0,arg1,(TM2-TM1)/1000.0); }

///////////////////////////////////////////////////////////////////////////////////////

/*
#include "stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <allegro.h>
#include <winalleg.h>
*/


//long NN[] = { 0, NZ , NY, NX };

////long NN[4];

#define SEEDS 20
//int     SeedVal;

int     SetEqualPower = 0;  
////float   StdSExp;//  = 1.0;
float   CmpSExp  = 1.0;
////float   StdTExp;//  = 1.0;
float   CmpTExp  = 1.0; 

////float  *FFT_Seed[SEEDS];         /* FFT of White Noise               */
////float  *FFT_Data;                /* Processing space for std or cmp  */





void fourn(float data[], long nn[], int ndim, int isign);

void randomize(void)
{
    time_t t;
    t = time(NULL);
    srand(t);
} 

void white_noise(float cmplx[],int len)
{
    int i;
    for(i=0; i<len; i++)
    {
        cmplx[i*2+0] = (float) rand();
        cmplx[i*2+1] = 0;
    }
}

void scale(float cmplx[],int len,float by)
{
    int i;
    for(i=0; i<len; i++)
    {
        cmplx[i*2+0] *= by;
        cmplx[i*2+1] *= by;
    }
}

void adjust_power(float cmplx[],int len)
{
    int i;
    float a,b,p,f;
    for(i=0; i<len; i++)
    {
        a = cmplx[i*2+0];
        b = cmplx[i*2+1];
        p = sqrt( a * a + b * b );
        f = atan2( b , a );
        p = 1;
        cmplx[i*2+0] = p*cos(f);
        cmplx[i*2+1] = p*sin(f);
    }
}

void filter(float cmplx[],float Exps,float Expt)
{
  int x,y,z;
  float fx,fy,fs,fz,ft;
  unsigned long i;
  long double m;

  i = 0;
  for (z=0; z<NZ; z++)
  {
      fz = (z<NZ/2)? z : NZ-z;
	  fz = fz * NZ_MULT;
      if (fz != 0) 
          ft = pow(fz,Expt);
      else 
          ft = 999999;
      for(y=0;y<NY;y++)
        {
          fy = (y<NY/2)? y : NY-y;
          for(x=0;x<NX;x++)
            {
              fx = (x<NX/2)? x : NX-x;
              fs = sqrt(fx*fx + fy*fy);

              if (fs==0)
                  m = 1;
              else
                  m = pow(fs,Exps);

              if (fz!=0)
                  m = m*ft;

              cmplx[i*2+0] = cmplx[i*2+0] / m;
              cmplx[i*2+1] = cmplx[i*2+1] / m;
              i++;
           }
        }
    }
}

void complex_to_float(float cmplx[],int len)
{
    int i;
    for (i=0; i<len; i++)
    {
        cmplx[i*2+0] = sqrt( cmplx[i*2+0]*cmplx[i*2+0]+cmplx[i*2+0]*cmplx[i*2+0]);
        cmplx[i*2+1] = 0;
    }
}

void InitData(void)
{
    FFT_Data = (float *)malloc(2 * NSIZE * sizeof(float));//new float[2*NSIZE];//(float *)calloc(2 * NSIZE, sizeof(float));
    srand(SeedVal);   
    white_noise(FFT_Data,NSIZE);
    fourn(FFT_Data-1, NN, 3, 1);
    scale(FFT_Data, NSIZE, 1.0 / NSIZE);    
    if (SetEqualPower) // set to false
        adjust_power(FFT_Data,NSIZE); 
    randomize();
}

void PrepareStimulus(float *fft_noise,float Exps,float Expt)
{
    filter(FFT_Data,Exps,Expt);
    fourn(FFT_Data-1,NN,3,-1); 
}



/*+***************************************************************
****
**    
****
***************************************************************+*/


#undef SWAP
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
/* Data[i*2 + 0] = ith re,
   Data[i*2 + 1] = ith im */
/* nn[] = dim for each */
/* ndim = for i=0; i<ndim;i++ => nn[i] */

void fourn(float data[], long nn[], int ndim, int isign)
{
    int idim;
    unsigned long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
    unsigned long ibit,k1,k2,n,nprev,nrem,ntot;
    float tempi,tempr;
    double theta,wi,wpi,wpr,wr,wtemp;

    for (ntot=1,idim=1;idim<=ndim;idim++)
        ntot *= nn[idim];
    nprev=1;
    for (idim=ndim;idim>=1;idim--) {
        n=nn[idim];
        nrem=ntot/(n*nprev);
        ip1=nprev << 1;
        ip2=ip1*n;
        ip3=ip2*nrem;
        i2rev=1;
        for (i2=1;i2<=ip2;i2+=ip1) {
            if (i2 < i2rev) {
                for (i1=i2;i1<=i2+ip1-2;i1+=2) {
                    for (i3=i1;i3<=ip3;i3+=ip2) {
                        i3rev=i2rev+i3-i2;
                        SWAP(data[i3],data[i3rev]);
                        SWAP(data[i3+1],data[i3rev+1]);
                    }
                }
            }
            ibit=ip2 >> 1;
            while (ibit >= ip1 && i2rev > ibit) {
                i2rev -= ibit;
                ibit >>= 1;
            }
            i2rev += ibit;
        }
        ifp1=ip1;
        while (ifp1 < ip2) {
            ifp2=ifp1 << 1;
            theta=isign*6.28318530717959/(ifp2/ip1);
            wtemp=sin(0.5*theta);
            wpr = -2.0*wtemp*wtemp;
            wpi=sin(theta);
            wr=1.0;
            wi=0.0;
            for (i3=1;i3<=ifp1;i3+=ip1) {
                for (i1=i3;i1<=i3+ip1-2;i1+=2) {
                    for (i2=i1;i2<=ip3;i2+=ifp2) {
                        k1=i2;
                        k2=k1+ifp1;
                        tempr=(float)wr*data[k2]-(float)wi*data[k2+1];
                        tempi=(float)wr*data[k2+1]+(float)wi*data[k2];

                        data[k2]=data[k1]-tempr;
                        data[k2+1]=data[k1+1]-tempi;

                        data[k1] += tempr;
                        data[k1+1] += tempi;
                    }
                }
                wr=(wtemp=wr)*wpr-wi*wpi+wr;
                wi=wi*wpr+wtemp*wpi+wi;
            }
            ifp1=ifp2;
        }
        nprev *= n;
    }
}

#undef SWAP
/* Fourn (C) Copr. 1986-92 Numerical Recipes Software %R&4. */
  
