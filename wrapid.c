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
	void Df(int plane,int ring, int sign, double **DD, fftw_plan pd) ;
	void calCR(double  *k_perp_Low, double *k_perp_High, double *k_perp_avg);
	void read_rings(int *rcountPtr, double **rdataPtr);
	void wrapid();
	void wrapid_fast();
	
	
	int r,i,j,k,l,ia,ja,ii,jj,r1,r2,r3,p1,p2,p3,rp,sign_2,sign_3,index, loop_id=0;
	int K1_bin, T_bin, MU_bin, NB=10, ell, emm,LL_k1,UL_k1,BW_k1;
	double Bispec_real, Num, tt, mu, Bispec_imag, Biana;
	int rank, size;
	
	
	double kk,dk_perp,time,start_time,dummy,dummy2;
	FILE *fp1, *fp2, *fp3,*fp13;
    char file1[128], file2[128], file3[128],file13[128], file_calcR[128],file_wrapid[128];
	
/*  Initializing Variables  */
        if(argc!=3){
        printf("Usage: %s NN(384) NRi(1)\n", argv[0]);
        return 1;}
	sscanf(argv[1],"%ld",&NN);
	sscanf(argv[2],"%d",&NRi);
	printf("\n Calculating BS Multipoles for Realization(NRi) = %d and Grids(NN) = %ld\n",NRi,NN);
	LL=.56;
	Nr=NN/3;
	Nmax=Nr;
	Nthreads=4;
	
	vol=pow((NN*LL),3);
	Cunits=2.*pi/(NN*LL);  
	
	ff=.1;
	beta=1.;
	


printf("\n Arrays are defined successfully\n");
double  *k_perp_Low, *k_perp_High, *k_perp_avg;
k_perp_Low=(double*)calloc(Nr,sizeof(double));
k_perp_High=(double*)calloc(Nr,sizeof(double));
k_perp_avg=(double*)calloc(Nr,sizeof(double));
dk_perp=(1.*Nmax)/(1.*Nr); 
  
  for(i=0;i<Nr;++i)
  {
   k_perp_Low[i]=i*dk_perp;
   k_perp_High[i]=(i+1.)*dk_perp;
   k_perp_avg[i]=(2./3.)*(pow(k_perp_High[i],3.)-pow(k_perp_Low[i],3.))/(pow(k_perp_High[i],2.)-pow(k_perp_Low[i],2.)); // bin average value
  }
k_perp_avg[0]=0.;

	RR= allocate_3d_double(Nr,NN,NN);

	for (r=0;r<Nr;++r)
	  for(i=0;i<NN;++i)
	    for(j=0;j<NN;++j)
	    {	
	      	ia=(i>NN/2) ? (i-NN) : i ;
		ja=(j>NN/2) ? (j-NN) : j ; 
		kk=sqrt(1.*(ia*ia+ja*ja));
		if(kk>=k_perp_Low[r] && kk<k_perp_High[r])
			  {RR[r][i][j]=1.;} 
		else  {RR[r][i][j]=0.;}
			
		}
	
/*------------ Initializing fft plan and openmp ------------*/
	fftw_init_threads();
	fftw_plan_with_nthreads(Nthreads);
	omp_set_num_threads(Nthreads);	

	

/*------------ Calculation of rings forming closed triangles on plane perpendicular to LoS ------------*/	
printf("\n 1. Calculating ring combinations\n");
start_time=omp_get_wtime();

snprintf(file_calcR,128*sizeof(char),"ResFiles_N%ld/data_files/closed_rings_%ld",NN,NN);
if (access(file_calcR, F_OK) == -1) {
       calCR(k_perp_Low,k_perp_High,k_perp_avg);} // stores possible combinations of rings and corresponding number of triangles
else { printf("Files already exist. Skipping calCR function.\n");}

printf("\n 2. Calculation of ring combinations completed successfully (in %f sec)\n",omp_get_wtime()-start_time);


  	
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
 
 
 
 
/*------------ Calculate and store IFT of data in each ring and plane ------------*/	
 	printf("\n 3. Calculating IFT of field in each ring and plane ....\n"); 
 	start_time=omp_get_wtime();	
 	
 	snprintf(file_wrapid,128*sizeof(char),"ResFiles_N%ld/data_files/dataPlaneRing_%ld_%ld_%ld_R%d",NN, NN,Nmax-4,Nr-4,NRi);
	if (access(file_wrapid, F_OK) == -1) {
		//wrapid(); // writes D for all rings and planes
		wrapid_fast();}
	else { printf("Files already exist. Skipping Wrapid function.\n");}
	printf("\n 4. Calculation done successfully (in %f sec)\n",omp_get_wtime()-start_time);
//#####################################################################################


}












