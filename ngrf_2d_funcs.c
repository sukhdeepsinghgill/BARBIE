

#include<stdlib.h>
#include<stdio.h>
#include<fftw3.h>
#include<math.h>
#include<unistd.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<omp.h>


//------ Global variables ---------//
extern long NN,N1,N2,N3;
extern int Nthreads;
extern double LL,ff; 
extern unsigned long int seed;
extern double Cx,Cy,Cz,tpibyL,pi,Area,vol;
extern double ***ro,**ro_2d;
extern fftw_plan p_ro, q_ro;
extern fftw_plan p_ro_2d, q_ro_2d;

/*-------------- Array Allocation functions -----------------*/
/*-------------- Array Allocation functions -----------------*/
double ***allocate_3d_double(long N1,long N2,long N3)   // dynamicL memory allocation to 3D array using pointers
{
  long ii,jj;
  long asize,index;
  double ***phia, *phi;

  phia=(double ***)calloc(N1,sizeof(double **));


  for(ii=0;ii<N1;++ii)
      phia[ii]=(double **)calloc(N2,sizeof(double *));

  asize = N1*N2;
  asize = asize*N3;

  if(!(phi = (double *) calloc(asize,sizeof(double))))
    {
      printf("error in allocate_3d_double.\n");
      exit(0);
    }

  for(ii=0;ii<N1;++ii)
    for(jj=0;jj<N2;++jj)
      {
	index = N2*N3;
	index = index*ii + N3*jj;
	phia[ii][jj]=phi+ index;   // it is address for contiguous (elements of the array are adjacent to one another in memory with no gaps between them)  memory allottment// adress increments with loop// phi is any random pointer to location, adding index to phi increment in location
      }						// it is allocting memory location to 0th element of [i][j] array (={1,2,3,4,...}, here 0th element is 1; see notebook
  return(phia);
}

double **allocate_2d_double(long N1,long N2)
{
  long ii;
  long asize;
  double **phia, *phi;

  phia=(double **)calloc(N1,sizeof(double *));
  
  asize = N1*N2;

  if(!(phi = (double *) calloc(asize,sizeof(double))))
    {
      fprintf(stderr,"Error in allocate_2d_double.\n");
      exit(0);
    }

  for(ii=0;ii<N1;++ii)
	phia[ii]=phi+ N2*ii;

  return(phia);
}



double ****allocate_4d_double(long Nr,long N1,long N2,long N3)   // dynamicL memory allocation to 3D array using pointers
{
  long ii,jj,ll;
  long asize,index;
  double ****phia, *phi;

  phia=(double ****)calloc(Nr,sizeof(double ***));


  for(ii=0;ii<Nr;++ii)
      phia[ii]=(double ***)calloc(N1,sizeof(double **));
      
	
  for(ii=0;ii<Nr;++ii)
	for(jj=0;jj<N1;++jj)
	  phia[ii][jj]=(double **)calloc(N2,sizeof(double *));

  asize = Nr*N1*N2*N3;

  if(!(phi = (double *) calloc(asize,sizeof(double))))
    {
      printf("error in allocate_4d_double.\n");
      exit(0);
    }

  for(ii=0;ii<Nr;++ii)
    for(jj=0;jj<N1;++jj)
    	for(ll=0;ll<N2;++ll)
      {
	index = N1*N2*N3*ii + N2*N3*jj + N3*ll;
	phia[ii][jj][ll]=phi+ index;   // it is address for contiguous (elements of the array are adjacent to one another in memory with no gaps between them)  memory allottment// adress increments with loop// phi is any random pointer to location, adding index to phi increment in location
      }						// it is allocting memory location to 0th element of [i][j] array (={1,2,3,4,...}, here 0th element is 1; see notebook
  return(phia);
}







//return input power spectrum
double P_T(double UU)// return input brightness temperature PS in K^2
{
  double ell,pu;  
  double Av=1.0, betav=-1.0; // amplitude and exponent of PS

  ell=2.0*M_PI*UU;// calculate ell
  pu=Av*pow(ell/1000.0,betav);      // 2pi/1000 fact
  return(pu);
}

double  p(double kk)
// return input brightness temperature PS in K^2
{
  double pu;  

  double Av=1.0, betav=-2.0; // amplitude and exponent of PS 
  pu=Av*pow(kk,betav);
  return(pu);
}
//_______________________________________________________//

double  Pk(long i, long j, long k)
// return input brightness temperature PS in K^2
{
  double pu,kx,ky,kz,kk;  

  kx=Cx*i;
  ky=Cy*j;
  kz=Cz*k;
  kk=sqrt(kx*kx+ ky*ky + kz*kz); 
  pu=p(kk);
  return(pu);
}
/*-------------------------------------------------------------*/





void power_spec(int Nbin, double *kmode, double *power,long *no)
{
  // local variables 
  
  long i, j, k, a, b, c;
  int d;
  fftw_complex *comp_ro;
  double fac1,fac2,fac3, fac, m, scale;
  long index,index1,index2,c1=0;
  
  comp_ro = (fftw_complex *)&(ro[0][0][0]);
  
  /*********** TO FIND POWER SPECTRUM OF RO. **************/
  
  fac1=1./(1.*N1*N1);
  fac2=1./(1.*N2*N2);
  fac3=1./(1.*N3*N3);
    
  for(i=0;i<Nbin;++i)
    {
      power[i]=0.0;
      kmode[i]=0.0;
      no[i]=0;
    }
  
  // do logarithmic bins
  
  scale=log10(0.5*N1)/Nbin;
  
  // ----------------------
  
  /**************** BINNING POWER SPECTRA **********************/
  
  /*-------------------------- half lines ------------------------- */
  
  for(i=1;i<=N1/2;i++)
    for(j=0;j<=N2/2;j=j+N2/2)
      for(k=0;k<=N3/2;k=k+N3/2)
	{
	  a=i;
	  b=j;
	  c=k;
	  
	  index = i*N2*(N3/2+1) + j*(N3/2+1) + k;
	  
	  m = sqrt(fac1*a*a + fac2*b*b + fac3*c*c);	      
	  
	  d=(int)floor(log10(m*N1)/scale);//logarithmic bins
	  
	  if(d>=0 && d<Nbin)
	    {
	      power[d]+= ((comp_ro[index][0]* comp_ro[index][0]) + (comp_ro[index][1]* comp_ro[index][1]))*1.;
	      kmode[d] += m*1.;
	      no[d]=no[d]+1;
	    }
	}	  
  
  /*----------------------- half planes -----------------------*/
  
  for(i=0;i<N1;i++)
    {
      a=(i>N1/2)? N1-i: i;
      index1 = i*N2*(N3/2+1) ;
      
      for(j=1;j<N2/2;j++) 
	{
	  b=j; 
	  index2 = index1 + j*(N3/2+1) ;
	  
	  for(k=0;k<=N3/2;k=k+N3/2)
	    {
	      c=k;
	      
	      index = index2 + k;
	      
	      m = sqrt(fac1*a*a + fac2*b*b + fac3*c*c);	      
	      
	      d=(int)floor(log10(m*N1)/scale);//logarithmic bins
	      
	      if(d>=0 && d<Nbin)
		{
		  power[d]+= ((comp_ro[index][0]* comp_ro[index][0]) + (comp_ro[index][1]* comp_ro[index][1]))*1.;
		  kmode[d] += m*1.;
		  no[d]=no[d]+1;
		}
	      
	    }	  
	}
    }
  
  /**************** half cube **********************/
  
  for(i=0;i<N1;i++)
    {
      a=(i>N1/2)? N1-i: i;
      index1 = i*N2*(N3/2+1) ;
      
      for(j=0;j<N2;j++)
	{
	  b=(j>N2/2)? N2-j: j;
	  index2 = index1 + j*(N3/2+1) ;
	  
	  for(k=1;k<N3/2;k++)
	    {
	      c=k;	  	      
	      
	      index = index2 + k;
	      
	      m = sqrt(fac1*a*a + fac2*b*b + fac3*c*c);	      
	      
	      /* m*(2 * pi/LL) is |k| */
	      /* m=1/2 corresponds to kmode[Nbin-1] i.e. Nyquits */
	      
	      
	      d=(int)floor(log10(m*N1)/scale);//logarithmic bins
	      
	      if (d>=0 && d<Nbin)
		{
		  power[d] += ((comp_ro[index][0]* comp_ro[index][0]) + (comp_ro[index][1]* comp_ro[index][1]))*1.;
		  kmode[d] += m*1.;
                no[d]=no[d]+1;
		}
	    } /* end k for loop */
	  
	}
    }
  
  
  
  for(i=0;i<Nbin;i++)
    {
      if (no[i]>0)
	{
	  power[i] = power[i]/(1.*no[i]*vol);
	  kmode[i]=tpibyL*kmode[i]/(1.*no[i]);
	}
    }
}
//=======================================================================//


void shapes(int nbin, double *uavg)
{
  int i,j,l; 
  double MU,TT, uu1;
  char file[128];
  FILE *fp1;
 /*----- Evaluating all possible unique shapes of closed triangles -----*/  
  for(i=0;i<nbin;++i)
  {
   uu1=Cx*uavg[i];
   
   snprintf(file, 128*sizeof(char),"shape_u%d_%.3lf",i,uu1);
   fp1=fopen(file,"w");
   
   for(j=0;j<=i;++j)
    {
     TT=uavg[j]/uavg[i];
     if(TT>=0.5 && TT<=1.)
     {
      for(l=0;l<=j;++l)
      {
       MU=0.5*(TT+1./TT-(uavg[l]/uavg[i])*(uavg[l]/uavg[j]));
       if(MU>=0.5/TT && MU<=1.) // this incorporates 2*mu*t>=1
          fprintf(fp1,"%e\t%e\n",TT,MU); // writing output
      }
     }
    }
   fclose(fp1);
  } 
}







//---- Given input power spectrum -------//
//---- fills modes for a Gaussian random field.......//
void Fill_Modes(unsigned long int seed) 
{  
  fftw_complex *A;
  A=(fftw_complex*)&(ro[0][0][0]);
  
  // local variables 
  long i,j,k,m,jj1,jj2,kk;
  long index1,index2;
  double amp,volfac;
  
  volfac=vol/2.;
  
  //for random number generator
  gsl_rng *r;
  double sig=1.;

  r= gsl_rng_alloc(gsl_rng_mt19937);
 gsl_rng_set (r, seed);
  
  //done

  /********************FILLING POINTS***********************/
  
  for(i=0;i<2;++i)
    for(j=0;j<2;++j)
      for(k=0;k<2;++k)
	{
	  if (i==0 && j==0 && k==0)
	    {
	      index1=0;
	      A[index1][0]=0.0;
	      A[index1][1]=0.0;
	    }
	  else
	    {
	      index1 = i*(N1/2)*N2*(N3/2+1) + j*(N2/2)*(N3/2+1) + k*(N3/2);
	      
              amp=sqrt(volfac*Pk(i*N1/2,j*N2/2,k*N3/2));
	      A[index1][0]=amp*gsl_ran_gaussian_ziggurat(r,sig);
	      A[index1][1]=0.0;
	    }
	}
  
  /********************FILLING LINES***********************/
  
  for(i=0;i<2;++i)
    for(j=0;j<2;++j)
      {
	for(k=1;k<N1/2;k++)
	  {
	    amp=sqrt(volfac*Pk(k,i*N2/2,j*N3/2));
	    
	    index1 =k*N2*(N3/2+1) +i*(N2/2)*(N3/2+1) + j*(N3/2);
	    
	    A[index1][0]=amp*gsl_ran_gaussian_ziggurat(r,sig);
	    A[index1][1]=amp*gsl_ran_gaussian_ziggurat(r,sig);
	    
	    index2 =(N1-k)*N2*(N3/2+1) +i*(N2/2)*(N3/2+1) + j*(N3/2);
	    
	    A[index2][0]=A[index1][0];
	    A[index2][1]=-A[index1][1];
	  }
	
	for(k=1;k<N2/2;k++)
	  {
	    amp=sqrt(volfac*Pk(i*N1/2,k,j*N3/2));
	    
	    index1 =i*(N1/2)*N2*(N3/2+1)+ k*(N3/2+1) + j*(N3/2);
	    index2 =i*(N1/2)*N2*(N3/2+1)+(N2-k)*(N3/2+1) + j*(N3/2);
	    
	    A[index1][0]=amp*gsl_ran_gaussian_ziggurat(r,sig);
	    A[index1][1]=amp*gsl_ran_gaussian_ziggurat(r,sig);
	    
	    A[index2][0]=A[index1][0];
	    A[index2][1]=-A[index1][1];
	  }
	
	for(k=1;k<N3/2;k++)
	  {
	    amp=sqrt(volfac*Pk(i*N1/2,j*N2/2,k));
	    
	    index1 =i*(N1/2)*N2*(N3/2+1)+ j*(N2/2)*(N3/2+1) + k;
	    
	    A[index1][0]=amp*gsl_ran_gaussian_ziggurat(r,sig);
	    A[index1][1]=amp*gsl_ran_gaussian_ziggurat(r,sig);
	  }
      }
  
  
  /**********************FILLING PLANES***********************/
  for(i=0;i<2;i++)
    for(m=0;m<2;m++)
      {
	for(j=1;j<N2/2;j++)
	  for(k=1;k<N3/2;k++)
	    {
	      amp=sqrt(volfac*Pk(i*N1/2,j,k));
	      
	      jj1=(1-m)*j + m * (N2-j);	    
	      index1 =i*(N1/2)*N2*(N3/2+1)+ jj1*(N3/2+1) + k;
	      
	      A[index1][0]=amp*gsl_ran_gaussian_ziggurat(r,sig);
	      A[index1][1]=amp*gsl_ran_gaussian_ziggurat(r,sig);	      
	    }
	
	for(j=1;j<N1/2;j++)
	  for(k=1;k<N3/2;k++)
	    {
	      amp=sqrt(volfac*Pk(j,i*N2/2,k));
	      
	      jj1=(1-m)*j + m * (N1-j);
	      index1 =jj1*N2*(N3/2+1)+ i*(N2/2)*(N3/2+1) + k;
	      
	      A[index1][0]=amp*gsl_ran_gaussian_ziggurat(r,sig);
	      A[index1][1]=amp*gsl_ran_gaussian_ziggurat(r,sig);
	    }	      
	
	for(j=1;j<N1/2;j++)
	  for(k=1;k<N2/2;k++)
	    {
	      amp=sqrt(volfac*Pk(j,k,i*N3/2));
	      
	      jj1=(1-m)*j + m * (N1-j);
	      kk=k;
	      index1 =jj1*N2*(N3/2+1)+ kk*(N3/2+1)+i*N3/2;
	      
	      jj2=m*j + (1-m) * (N1-j);
	      kk=N2-k;
	      index2 =jj2*N2*(N3/2+1)+ kk*(N3/2+1)+i*N3/2;
	      
	      A[index1][0]=amp*gsl_ran_gaussian_ziggurat(r,sig);
	      A[index1][1]=amp*gsl_ran_gaussian_ziggurat(r,sig);
	      
	      A[index2][0]=A[index1][0];
	      A[index2][1]=-A[index1][1];
	    }
      }
  /***************************FILLING CUBES***********************/ 
  
  for(i=1;i<N1/2;i++)
    for(j=1;j<N2/2;j++)
      for(k=1;k<N3/2;k++)
	{
	  amp=sqrt(volfac*Pk(i,j,k));
	  
	  index1 =i*N2*(N3/2+1) + j*(N3/2+1) + k;
	  A[index1][0]=amp*gsl_ran_gaussian_ziggurat(r,sig);
	  A[index1][1]=amp*gsl_ran_gaussian_ziggurat(r,sig);
	  
	  index1 =(N1-i)*N2*(N3/2+1) + j*(N3/2+1) + k;
	  A[index1][0]=amp*gsl_ran_gaussian_ziggurat(r,sig);
	  A[index1][1]=amp*gsl_ran_gaussian_ziggurat(r,sig);
	  
	  index1 =i*N2*(N3/2+1) + (N2-j)*(N3/2+1) + k;
	  A[index1][0]=amp*gsl_ran_gaussian_ziggurat(r,sig);
	  A[index1][1]=amp*gsl_ran_gaussian_ziggurat(r,sig);
	  
	  index1 = (N1-i)*N2*(N3/2+1) + (N2-j)*(N3/2+1) + k;
	  A[index1][0]=amp*gsl_ran_gaussian_ziggurat(r,sig);
	  A[index1][1]=amp*gsl_ran_gaussian_ziggurat(r,sig);
	}

  // finished filling Fourier modes
 gsl_rng_free(r);
}





/*------- Changes Fourier data to image --------*/
void Make_Image() 
{ 
  // local variables
  int i,j,k,index;
  
  // for filling Fourier components 
  fftw_complex *A;
  A=(fftw_complex*)&(ro[0][0][0]);
  
  for(i=0;i<N1;i++)
    for(j=0;j<N2;j++)
      for(k=0;k<=N3/2;k++)
	 {
	  index=i*N2*(N3/2+1) + j*(N3/2+1) + k;
	  A[index][0]=A[index][0]/(vol);
	  A[index][1]=A[index][1]/(vol);
        }
  
  fftw_execute(q_ro); // perform inverse FT, ro now contains delta(x)
}

/*-------- Getting Fourier data from image ---------*/
void Get_Modes()
// Image --> Fourier Modes
{
  // local variables
  int i,j,k;
  double Lcube=pow(LL,3.);
  
  for(i=0;i<N1;i++)
    for(j=0;j<N2;j++)
      for(k=0;k<N3;k++)
	  ro[i][j][k]=ro[i][j][k]*Lcube;

      
      fftw_execute(p_ro);  // ro now contains Delta(k)

}
/*----------------------------------------------------------*/

/*------- Introducing non-Gaussianity to the image ---------*/
double Non_Gauss()
// Given Gaussian random image, this introduces non-Gaussianity 
// d=d + fNL*(d^2-var) 
{ 
  double var=0.; // var--> variance of Gaussian field
  int i,j,k;
  
 //---- calculate variance ----// 
  for(i=0;i<N1;++i)
    for(j=0;j<N2;++j)
      for(k=0;k<N3;++k)
	  var+=pow(ro[i][j][k],2.);

  var/=(1.*N1*N2*N3);
 //---- done variance ----//
   
 //--- non-Gaussianity is being introduced ---//
 for(i=0;i<N1;++i)
    for(j=0;j<N2;++j)
      for(k=0;k<N3;++k)
        ro[i][j][k]+=ff*(pow(ro[i][j][k],2.)-var);
 //--- non-Gaussianity introduced ---//
        
 return(pow(var,0.5)); // return sigma in the main program
}









void Fill_Modes_2D(unsigned long int seed)
{  
  // for filling Fourier components 
  fftw_complex *in;
  in=(fftw_complex*)&ro_2d[0][0];

  // local variables 
  double Length, uu, amp, fact;
  int i,j,ia,index,index1,ydim;
  
  Length=NN*LL;
  fact=Length/sqrt(2.0); // factor for amplitude
  ydim=(N2/2+1); 
  
  //for random number generator
  gsl_rng *r;
  double sigma=1.;

  r= gsl_rng_alloc(gsl_rng_cmrg);
  gsl_rng_set (r, seed);
  
  /*------ Filling Fourier Components --------*/ 

  //along axis (j-0 and j=N/2)

  for(j=0;j<ydim;j=j+N2/2)
    for(i=1;i<N1/2;++i)
      {
	// along + x 
	uu=sqrt(1.*(i*i+j*j))/Length;
	amp=fact*sqrt(P_T(uu));
	index=i*ydim+j;
	
	in[index][0]=pow(-1.,i+j)*amp*gsl_ran_gaussian(r,sigma);
	in[index][1]=pow(-1.,i+j)*amp*gsl_ran_gaussian(r,sigma);
	// along -x 
	index1=(N1-i)*ydim+j;
	
	in[index1][0]=in[index][0];
	in[index1][1]=-in[index][1];
      }
  // upper half plane excluding x axis
  for(i=0;i<N1;++i)
    for(j=1;j<N2/2;++j)
      {
	ia= (i>N1/2) ? (N1-i) : i ;
	uu=sqrt(1.*(ia*ia+j*j))/Length;
	amp=fact*sqrt(P_T(uu));
	index=i*ydim+j;
	
	in[index][0]=pow(-1.,i+j)*amp*gsl_ran_gaussian(r,sigma);
	in[index][1]=pow(-1.,i+j)*amp*gsl_ran_gaussian(r,sigma);
      }
      
  // remianing 4 points
  for(i=0;i<2;++i)
    for(j=0;j<2;++j)
      {
	if(i+j==0) 
	  {
	    in[0][0]=0.0;
	    in[0][1]=0.0;
	  }
	else
	  {
	    uu=(NN/2.)*sqrt(1.*(i*i+j*j))/Length;
	    amp=fact*sqrt(P_T(uu));
	    index=i*(N1/2)*ydim+j*(N2/2);
	    
	    in[index][0]=pow(-1.,(i*N1/2+j*N2/2))*amp*gsl_ran_gaussian(r,sigma);
	    in[index][1]=0.0;
	  }
      }
 // finished filling Fourier components
}

// Fourier Modes -> Image
void Make_Image_2D()
{ 
  // local variables
  int i,j,index;

  // for filling Fourier components 
  fftw_complex *in;
  in=(fftw_complex*)&ro_2d[0][0];
  
  for(i=0;i<N1;++i)
    for(j=0;j<=N2/2;++j)
      {
	index=i*(N2/2+1)+j;
	in[index][0]/=Area;
	in[index][1]/=Area;
      }
  fftw_execute(q_ro_2d); // perform inverse FT, ro now contains delta(x)
}
// Finished

// Image --> Fourier Modes
void Get_Modes_2D()
{
  // local variables
  int i,j,index;
  double LLsq=pow(LL,3.0);

  for(i=0;i<N1;++i)
    for(j=0;j<N2;++j)
	ro_2d[i][j]*=LLsq;

  fftw_execute(p_ro_2d);
}
// Finished

// Given Gaussian random image, this introduces non-Gaussianity 
// d=d + ff*(d*2-var)
double Non_Gauss_2D()
{ 
  double var=0; // var --> variance of Gaussian field
  int i,j,index;
  
  // Computing variance 
  for(i=0;i<N1;++i)
    for(j=0;j<N2;++j)
	var+=pow(ro_2d[i][j],2.);
      
  var/=(1.*N1*N2);

 // Introducing non-Gaussianity 
 for(i=0;i<NN;++i)
    for(j=0;j<NN;++j)
	ro_2d[i][j]+=ff*(pow(ro_2d[i][j],2.)-var);//*pow(var,-0.5); // non-Gaussianity introduced

 return(pow(var,0.5)); // return sigma in the main program to compute fNL/sigma
}













/*------ Introduces linear-order redshift space distortion in data -------*/
void RSD(double beta)
{ // beta is the RSD parameter

 //----- Initializing variables -----//
 double myu,ia,ja,kk;
 int i,j,k,index;
 fftw_complex *A;

 A=(fftw_complex *)&ro[0][0][0];
 
 //---- modifying Delta(k)=(1+beta*mu^2)*Delta(k) ----//
 #pragma omp parallel for num_threads(Nthreads) collapse(3) default(shared) private(myu,index,ia,ja)
 for(i=0;i<N1;++i)
   for(j=0;j<N2;++j)
     for(k=0;k<=N3/2;++k)  
       {
       	
    ia=(i>N1/2) ? (N1-i) : i ;   // condition ? true : false
  	ja=(j>N2/2) ? (N2-j) : j ;
  	index= i*N2*(N3/2+1) + j*(N3/2+1) + k;  // fourier transform harmitian half sphere // due to complex array handling
  	
	kk=sqrt(1.*(ia*ia+ja*ja+k*k));
	 if(i+j+k!=0)
	 {
	  myu=k/kk; // LoS component of k vector    why k/sqrt(i2+j2+k2) is LOS???????????? choose k is los
	  index=i*N2*(N3/2+1)+ j*(N3/2+1)+ k;
	  A[index][0]*=(1.+beta*myu*myu);
	  A[index][1]*=(1.+beta*myu*myu);
        }
      }
 //------ done ------//
}















