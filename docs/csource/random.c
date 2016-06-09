/* $Id$ */
/*
 * random.c
 *
 * written by Sebastian Voecking <seb.voeck@uni-muenster.de>
 *
 * For details see random.h
 */

#include "random.h"

#include <stdlib.h>

static int method = RANDOM_STDLIB;

static int rand_seed=-1; // for random calculation
static double get_random_1_cw();
static void subrn(double* u, int len);
static double random_james();

void random_set_method(RandomMethode meth)
{
  method = meth;
}

double random_get()
{
  switch(method) {
    case RANDOM_STDLIB:
      return ((double)rand())/RAND_MAX;
    case RANDOM_CW:
      return get_random_1_cw();
    case RANDOM_JAMES:
      return random_james();
    default:
      return 0;
  }
}

void random_seed(int seed)
{
  switch(method) {
    case RANDOM_STDLIB:
      srand(seed);
      break;
    case RANDOM_CW:
      rand_seed = seed;
      break;
    case RANDOM_JAMES:
      break;
  }
}

double get_random_1_cw()
     // function gives a random value between 0 and 1.
     // NO imput value needed.
     // numerical rec. version provided by Ch. Weinheimer
{
#define IA 16807
#define IM 2147483647
#define AM (1./IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define RNMX (1.-1.2e-7)
  
  int j, k;
  static int iy=0, iv[NTAB];
  double temp;
  
  if (rand_seed<=0 || !iy){
    if (-rand_seed<1) rand_seed=1;
    else rand_seed=-rand_seed;
    for (j=NTAB+7; j>=0; j--){
      k=rand_seed/IQ;
      rand_seed=IA*(rand_seed-k*IQ)-IR*k;
      if (rand_seed<0) rand_seed+=IM;
      if (j<NTAB) iv[j]=rand_seed;
    }
    iy=iv[0];
  }
  k=rand_seed/IQ;
  rand_seed=IA*(rand_seed-k*IQ)-IR*k;
  if (rand_seed<0) rand_seed+=IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j]=rand_seed;
  if ((temp=AM*iy)>RNMX) return RNMX;
  else return temp;
  
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef RNMX
}

void subrn(double *u,int len)
{
// This subroutine computes random numbers u[1],...,u[len]
// in the (0,1) interval. It uses the 0<IJKLRANDOM<900000000
// integer as initialization seed.
//  In the calling program the dimension
// of the u[] vector should be larger than len (the u[0] value is
// not used).
// For each IJKLRANDOM
// numbers the program computes completely independent random number
// sequences (see: F. James, Comp. Phys. Comm. 60 (1990) 329, sec. 3.3).
  //
  // remark by T. Thuemmler:
  // same random numbers appear each time one restarts the whole program
  //
  static long IJKLRANDOM=100;
  static int iff=0;
  static long ijkl,ij,kl,i,j,k,l,ii,jj,m,i97,j97,ivec;
  static float s,t,uu[98],c,cd,cm,uni;
  if(iff==0)
  { ijkl=IJKLRANDOM;
    if(ijkl<1 || ijkl>=900000000) ijkl=1;
    ij=ijkl/30082;
    kl=ijkl-30082*ij;
    i=((ij/177)%177)+2;
    j=(ij%177)+2;
    k=((kl/169)%178)+1;
    l=kl%169;
    for(ii=1;ii<=97;ii++)
    { s=0; t=0.5;
      for(jj=1;jj<=24;jj++)
      { m=(((i*j)%179)*k)%179;
        i=j; j=k; k=m;
        l=(53*l+1)%169;
        if((l*m)%64 >= 32) s=s+t;
        t=0.5*t;
      }
      uu[ii]=s;
    }
    c=362436./16777216.;
    cd=7654321./16777216.;
    cm=16777213./16777216.;
    i97=97;
    j97=33;
    iff=1;
  }
  for(ivec=1;ivec<=len;ivec++)
  { uni=uu[i97]-uu[j97];
    if(uni<0.) uni=uni+1.;
    uu[i97]=uni;
    i97=i97-1;
    if(i97==0) i97=97;
    j97=j97-1;
    if(j97==0) j97=97;
    c=c-cd;
    if(c<0.) c=c+cm;
    uni=uni-c;
    if(uni<0.) uni=uni+1.;
    if(uni==0.)
    { uni=uu[j97]*0.59604644775391e-07;
      if(uni==0.) uni=0.35527136788005e-14;
    }
    u[ivec]=uni;
  }

  //  cout << endl<<  "random: " << u[1] << endl << flush;

  return;
}

double random_james()
{
// This function computes 1 random number in the (0,1) interval,
// using the subrn subroutine.
  double u[2];
  subrn(u,1);
  return u[1];
}

///////////////////////////////////////////////////////////



