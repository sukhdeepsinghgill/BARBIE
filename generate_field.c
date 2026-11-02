#include<stdlib.h>
#include<stdio.h>
#include<fftw3.h>
#include<math.h>
#include<unistd.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<omp.h>


long NN,N1,N2,N3;
int Nthreads, F_Task;
double LL,ff,beta;
unsigned long int seed;
double Cx,Cy,Cz,tpibyL,Area,vol;
double ***ro,**ro_2d;
double pi=M_PI; 
fftw_plan p_ro, q_ro;
fftw_plan p_ro_2d, q_ro_2d;


int main(int argc, char *argv[])
{
	/*  Declaring Variables  */
	
	double ***allocate_3d_double(long N1,long N2,long N3);
	double **allocate_2d_double(long N1,long N2);
	void Fill_Modes(unsigned long int seed);
	
	int r,i,j,k,index,F_RSD,NRi;
	

/*  Initializing Variables  */

        if(argc!=3){
        printf("Usage: %s <Number of Grids> <Realization>\n", argv[0]);
        return 1;}
	sscanf(argv[1],"%ld",&NN);
	sscanf(argv[2],"%d",&NRi);
	
	printf("\n Generating field Realization(NRi) = %d for Grids(NN) = %ld\n\n",NRi,NN);

	
	LL=0.56;
	Nthreads=4;
	
	seed=100020+146*NRi;
	ff=0.5;
	beta=1.;
	F_RSD=1;
	
	vol=pow((NN*LL),3);
	N1=NN;N2=NN;N3=NN; 
	tpibyL=2.*pi/LL;
  	Cx=2.*pi/(N1*LL);  Cy=2.*pi/(N2*LL);   Cz=2.*pi/(N3*LL); 

	Nthreads=4;
	
	  
  

  
 /*------------ Initializing fft plan and openmp ------------*/
  fftw_init_threads();
  fftw_plan_with_nthreads(Nthreads);
  omp_set_num_threads(Nthreads);

  //memory allocation for density field
  ro= allocate_3d_double(N1,N2,N3+2);
//  ro_2d= allocate_2d_double(NN,NN*2);

  
  double sigma;
  void Make_Image();
  void Get_Modes();
  double Non_Gauss();
  void RSD(double beta);
  p_ro = fftw_plan_dft_r2c_3d(N1, N2, N3, &(ro[0][0][0]), (fftw_complex*)&(ro[0][0][0]), FFTW_ESTIMATE);  
  q_ro = fftw_plan_dft_c2r_3d(N1, N2, N3, (fftw_complex*)&(ro[0][0][0]), &(ro[0][0][0]), FFTW_ESTIMATE);
 
  Fill_Modes(seed);
  Make_Image();
  sigma = Non_Gauss();
  fprintf(stdout,"sigma = %e\n",sigma); 
  //ff/=sigma; 
  ///fprintf(stdout,"fNL/sigma = %e\n",ff);
  Get_Modes();
  RSD(beta);

  fftw_complex *comp_ro; 
   comp_ro = (fftw_complex *)&(ro[0][0][0]);

//  fftw_complex *comp_ro_2d;
//    comp_ro_2d = (fftw_complex *)&(ro_2d[0][0]);
    
    FILE *pp;
    char fname[128];
    snprintf(fname,128*sizeof(char),"ResFiles_N%ld/data_files/field_N%ld_R%d.txt",NN,NN,NRi);
    pp=fopen(fname,"w");
    
    for(i=0;i<N1;++i)
      for(j=0;j<N2;++j)
        for(k=0;k<N3+2;++k){
			fwrite(&ro[i][j][k],sizeof(double),1,pp); // we can also write: fwrite(&comp_ro[index][0],sizeof(double),1,pp); fwrite(&comp_ro[index][1],sizeof(double),1,pp); where index=i*N2*(N3/2+1) + j*(N3/2+1)  + k; and loop k=0;k<N3/2
       }
        
    fclose(pp);
    
//    for(i=0;i<N1;++i){
//      for(j=0;j<N2;++j){
//        for(k=0;k<N3;++k){
//        	index=i*N2*(N3/2+1) + j*(N3/2+1)  + k;
//         	//printf("%f,%f,",comp_ro[index][0],comp_ro[index][1]);
//			printf("%f,",ro[i][j][k]);
//		}printf("\n");
//	}printf("\n");
//}
//    
//    
    
    
    
    
}
