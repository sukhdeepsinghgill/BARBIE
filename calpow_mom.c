// This version is modification of zBS and written after zbS_v2. Here implimentation of loops is changed to make it even faster taan zBS (with FBI_fast) i.e implementation is as zBS but Loop order is as zBS_v2. This is slighly faster than zBS (165 to 155 sec etc)

#include<stdlib.h>
#include<stdio.h>
#include<fftw3.h>
#include<math.h>
#include<unistd.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<omp.h>

//------ Global variables ---------//
long NN,Nmax,Nr;
double LL,Cunits,vol;
int Nthreads, NRi;
double ***ro,**ro_2d, ***RR, **D1, **D2, **D3;
double pi=M_PI; 
fftw_plan pd1, pd2, pd3;



int main(int argc, char *argv[])
{
/*  Declaring Variables  */

	
	
	double dK, K1, T, MU, swap, px, pz, Ylm, norm_fact,mu_1,mu_2,mu_3,ff,beta;
	double ***allocate_3d_double(long N1,long N2,long N3);
	double **allocate_2d_double(long N1,long N2);
	double  input_ps(double kk);
	void calpow_mom(int Nbin,double* power, double* kmode, double* power_P2,double* power_P4, double* power_P6,double *no);
	
	
	int r,i,j,k,l,ia,ja,ii,jj,r1,r2,r3,p1,p2,p3,rp,sign_2,sign_3,index, loop_id=0;
	int K1_bin, T_bin, MU_bin, ell, emm,LL_k1,UL_k1,BW_k1;

	
	double kk,dk_perp,time,start_time,dummy,dummy2;
	FILE *fp1, *fp2, *fp3,*fp13;
    char file1[128], file2[128], file3[128];
	
/*  Initializing Variables  */
        if(argc!=3){
        printf("Usage: %s NN(384) NRi(1) UL_k1(Nmax) LL_k1(Nmax-10) BW_k1(1)\n", argv[0]);
        return 1;}
	sscanf(argv[1],"%ld",&NN);
	sscanf(argv[2],"%d",&NRi);
	printf("\n Calculating PS Multipoles for Realization(NRi) = %d and Grids(NN) = %ld\n",NRi,NN);
	LL=.56;
	Nmax=NN/3;
	Nr=Nmax/4;
	Nthreads=4;
	
	vol=pow((NN*LL),3);
	Cunits=2.*pi/(NN*LL);  
	
	ff=0.5;
	beta=1.;
	

	
/*------------ Initializing fft plan and openmp ------------*/
	fftw_init_threads();
	fftw_plan_with_nthreads(Nthreads);
	omp_set_num_threads(Nthreads);	

  	
/*------------ Read Field ------------*/
  	ro= allocate_3d_double(NN,NN,NN+2);
   
  	snprintf(file1,128*sizeof(char),"ResFiles_N%ld/data_files/field_N%ld_R%d.txt",NN,NN,NRi);
   	fp1=fopen(file1,"r");
    	for(i=0;i<NN;++i)
      	  for(j=0;j<NN;++j)
            for(k=0;k<NN+2;++k){
		fread(&ro[i][j][k],sizeof(double),1,fp1); }
        
    	fclose(fp1);	
  	
	fftw_complex *comp_ro; 
	comp_ro = (fftw_complex *)&(ro[0][0][0]);
 
 


 /*-------------- power spectrum computation ---------------*/ 
  int NB=10;
  double *num;
  double *kmean,*pkmean, *pk2mean, *pk4mean,*pk6mean;
  double *powerk;
 
  num=(double*)calloc(NB,sizeof(long));
  kmean=(double*)calloc(NB,sizeof(double));
  pkmean=(double*)calloc(NB,sizeof(double));
  powerk=(double*)calloc(NB,sizeof(double));
  pk2mean=(double*)calloc(NB,sizeof(double));
  pk4mean=(double*)calloc(NB,sizeof(double));
  pk6mean=(double*)calloc(NB,sizeof(double));
  
  
  start_time=omp_get_wtime();

 
  calpow_mom(NB,pkmean, kmean, pk2mean,pk4mean,pk6mean, num);


  /*----------- Write power spectrum in data file -----------*/
  snprintf(file3, 128*sizeof(char),"ResFiles_N%ld/psoutN%ld_R%d.dat",NN,NN,NRi);
  fp3=fopen(file3,"w");

  for(i=0;i<NB;++i)
    if(num[i]>0)
      {
  	fprintf(fp3,"%e\t%e\t%e\t%e\t%e\t%e\n",kmean[i],pkmean[i],pk2mean[i],pk4mean[i],pk6mean[i],num[i]);
      }
  fclose(fp3);
  
  //fprintf(stdout,"Power spectrum time = %e sec\n",omp_get_wtime()-start_time);

}


  

