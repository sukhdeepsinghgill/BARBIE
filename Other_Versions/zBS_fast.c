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
	void FBI(int p1,int r1,int p2,int r2,int p3,int r3,int sign_2,int sign_3,double *Bispec_real, double *Bispec_imag);
	void FBI_fast(int p1,int r1,int p2,int r2,int p3,int r3,int sign_2,int sign_3,double *Bispec_real, double *Bispec_imag);
	void cal_coeff(double px, double pz, double *bs, double *norm);
	void sort(double K[3][2]);
	double Expo(int m, double px, double pz);
	double gsl_sf_legendre_sphPlm(int l, int m, double x);
	double gsl_sf_legendre_Plm (int l, int m, double x);
	
	int r,i,j,k,l,ia,ja,ii,jj,r1,r2,r3,p1,p2,p3,rp,sign_2,sign_3,index, loop_id=0;
	int K1_bin, T_bin, MU_bin, NB=10, ell, emm,LL_k1,UL_k1,BW_k1;
	double Bispec_real, Num, tt, mu, Bispec_imag, Biana;
	
	double kk,dk_perp,time,start_time,dummy,dummy2;
	FILE *fp1, *fp2, *fp3,*fp13;
    char file1[128], file2[128], file3[128],file13[128], file_calcR[128],file_wrapid[128];
	
/*  Initializing Variables  */
        if(argc!=6){
        printf("Usage: %s NN(384) NRi(1) UL_k1(Nmax) LL_k1(Nmax-10) BW_k1(1)\n", argv[0]);
        return 1;}
	sscanf(argv[1],"%ld",&NN);
	sscanf(argv[2],"%d",&NRi);
	sscanf(argv[3],"%d",&UL_k1); //Upper limit on K1
	sscanf(argv[4],"%d",&LL_k1); // Lower limit on K1
	sscanf(argv[5],"%d",&BW_k1); // Bin width
	printf("\n Calculating BS Multipoles for Realization(NRi) = %d and Grids(NN) = %ld\n",NRi,NN);
	LL=.56;
	Nthreads=4;


	Nmax=86; //150
	Nr=32;   //40
	dk_perp=1.; // i'll keep changing dkperp in loop & make sure dperplow[0]>=1;dperpHigh[Nr-1]<=Nmax
	
	/*Nmax=NN/3;
	Nr=40;//Nmax/4;
	dk_perp=(1.*Nmax)/(1.*Nr);*/ //remove 1 in kperplow&kperphigh and set kperplow[0]=1//or add 1 in kperplow&kperphigh & dk_perp=(Nmax-1)/Nr
	
	/*dk_perp = 4.;
	Nmax = (NN/3)-(dk_perp-1);
	Nr = (Nmax-1)/dk_perp;*/  // add 1 in kperplow and kperphigh

	/*Nmax = NN/3;
	Nr   = 20;
	dk_perp=(1.*log10(Nmax))/(1.*Nr);*/  // for log binning
	
	vol=pow((NN*LL),3);
	Cunits=2.*pi/(NN*LL);  
	
	ff=0.5;
	beta=1.;
	


int n_mtbin=6, npoles=16, n_k1bin= (int)((UL_k1-LL_k1)/(1.*BW_k1)) ;//calc will be saved for last 10 bins of width 2 grids

printf("\n All parameters are with me including mu-t bins=%d, k1 bins=%d  !!\n",n_mtbin,n_k1bin);
double bin_BS_real[n_k1bin][n_mtbin][n_mtbin][npoles],bin_normalization[n_k1bin][n_mtbin][n_mtbin][npoles], bs[npoles], norm[npoles],bin_Biana[n_k1bin][n_mtbin][n_mtbin][npoles];//,bin_BS_imag[Nmax][n_mtbin][n_mtbin][npoles];  //bin_tri[Nmax][n_mtbin][n_mtbin]
printf("\n Arrays are defined successfully\n");
double  *k_perp_Low, *k_perp_High, *k_perp_avg;
k_perp_Low=(double*)calloc(Nr,sizeof(double));
k_perp_High=(double*)calloc(Nr,sizeof(double));
k_perp_avg=(double*)calloc(Nr,sizeof(double));

 /*
  for(i=0;i<Nr;++i)
  {
   k_perp_Low[i]=i*dk_perp ;
   k_perp_High[i]=(i+1.)*dk_perp ;
   //k_perp_avg[i]=(2./3.)*(pow(k_perp_High[i],3.)-pow(k_perp_Low[i],3.))/(pow(k_perp_High[i],2.)-pow(k_perp_Low[i],2.)); // bin average value
   k_perp_avg[i]=(k_perp_High[i]+k_perp_Low[i])/2.;
  }
*/
//k_perp_Low[0] = 1.;
//k_perp_avg[0]=0.;   




/*-------- varying dkperp in binning ----------*/

  int currentNumber = 0;
 /* for ( i = 0; i < Nr; i += 10) {
    for ( j = i; j < i + 10; j++) {
        k_perp_Low[j] = currentNumber;
	currentNumber += (i / 10) + 1;}}*/
    for (i = 0; i < Nr; i += 10) {
    int step_size = 1 << (i / 10);  
    for (j = i; j < i + 10; j++) {
        k_perp_Low[j] = currentNumber;
        currentNumber += step_size;}}  

for(i=1;i<Nr;++i)
  {
    k_perp_Low[i] -= .1;  // otherwise mu=1, stretched bs=nan; could do +=.1 but then bs for k=0 mode is difficult to handle.so keep 0th ring from 0-.9
    k_perp_High[i-1]=k_perp_Low[i];
  }
  k_perp_High[Nr-1]=1.*Nmax;
  
  for(i=0;i<Nr;++i)
  { k_perp_avg[i]=(k_perp_High[i]+k_perp_Low[i])/2.;

/*if (i<10)
   {k_perp_avg[i]=(k_perp_High[i]+k_perp_Low[i])/2.;}
else
   {k_perp_avg[i]=(2./3.)*(pow(k_perp_High[i],3.)-pow(k_perp_Low[i],3.))/(pow(k_perp_High[i],2.)-pow(k_perp_Low[i],2.)); }// bin average value */
   printf("\t%f\t%f\t%f\n",k_perp_Low[i],k_perp_High[i],k_perp_avg[i]) ;
  }

/*-------- end ------------*/



/*-------- log binning in dkperp ----------*/
/*  for(i=0;i<Nr;++i)
  {
   k_perp_Low[i]=pow(10,i*dk_perp) ;
   k_perp_High[i]=pow(10,(i+1.)*dk_perp) ;
   k_perp_avg[i]=(2./3.)*(pow(k_perp_High[i],3.)-pow(k_perp_Low[i],3.))/(pow(k_perp_High[i],2.)-pow(k_perp_Low[i],2.)); // bin average value
  }
*/




//	double K_min=k_perp_avg[0], K_max=sqrt(k_perp_avg[Nr-1]*k_perp_avg[Nr-1] + NN*NN/9);
//	dK=.05*log10((0.6*NN/k_perp_avg[0]));	
//	dK=.1*log10(K_max/K_min);	
/*  Arrays Assignment  */
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


/*------------ Reading data of rings forming closed triangles on plane perpendicular to LoS ------------*/	
int rcount;
double *rdata;
read_rings(&rcount, &rdata);  // read data generated by CalCR
/*printf("%d\n",rcount);
for (int i = 0; i < rcount; i++) {
 for (int j = 0; j < 4; j++) {
    printf("%lf\t ", rdata [i * 4 + j]);}
 printf("\n");}*/

/////////////////////////////////////////////////////////////////////////////
	
	/*D1=allocate_2d_double(NN,NN*2); 	
	D2=allocate_2d_double(NN,NN*2); 
	D3=allocate_2d_double(NN,NN*2);*/ 
	


  	
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

//#####################################################################################



/////////////////////////////////////////////////////////////////////////
for (i=0;i<n_k1bin;++i){
 for (j=0;j<n_mtbin;++j){
  for (l=0;l<n_mtbin;++l){
   for (jj=0;jj<npoles;++jj){
	bin_BS_real[i][j][l][jj]=0.;
	//bin_BS_imag[i][j][l][jj]=0.;
	bin_normalization[i][j][l][jj]=0.;} 
}}}
	

/*char filek1[128];
snprintf(filek1, 128*sizeof(char),"ResFiles_N%ld/all_k1_%ld",NN,NN);
FILE *fk1=fopen(filek1,"wb");*/

start_time=omp_get_wtime();
snprintf(file2, 128*sizeof(char),"ResFiles_N%ld/time_%ld_R%d.txt",NN,NN,NRi);	
snprintf(file1, 128*sizeof(char),"ResFiles_N%ld/FBiANDbiana_%ld_refined.txt",NN,NN);		
fp2 = fopen(file2, "w"); // Open the file in write mode to create a fresh file
fclose(fp2); 	
loop_id=0;
time=omp_get_wtime();
char drfilename1[128];
char drfilename2[128];
char drfilename3[128];
FILE *drfile2,*drfile3, *drfile1;
double buffer1[2*NN*NN];
double buffer2[2*NN*NN];
double buffer3[2*NN*NN];

int p2_min,p2_max,p3_min,p3_max;



/*
for(rp=0;rp<rcount;++rp)
{ 
  r1=(int)rdata [rp * 4 + 0] ; //rdata[rp][0];
  r2=(int)rdata [rp * 4 + 1] ; // rdata[rp][1];
  r3=(int)rdata [rp * 4 + 2] ; //rdata[rp][2];
  //Num=rdata [rp * 4 + 3] ; // rdata[rp][3];
   	 
for(p1=0;p1<Nmax;++p1)
{   */ 	

for(p1=0;p1<Nmax;++p1)
{ 
for(rp=0;rp<rcount;++rp)
{ 
  r1=(int)rdata [rp * 4 + 0] ; //rdata[rp][0];
  r2=(int)rdata [rp * 4 + 1] ; // rdata[rp][1];
  r3=(int)rdata [rp * 4 + 2] ; //rdata[rp][2];
  //Num=rdata [rp * 4 + 3] ; // rdata[rp][3];

	
if (r1==r2 && r2 ==r3){
		p2_min = 0  ;
		p2_max = p1  ;}
	else if (r1==r2 && r2 != r3){
		p2_min = -p1  ;
		p2_max = p1  ;}
	else if (r1 != r2 && r2 ==r3){
		p2_min =  (int)(-p1/2) ;
		p2_max =  Nmax ;}
	else{
		p2_min = (p1 == 0) ? 0 : -Nmax + 1;  
		p2_max = Nmax  ;}
	

	/*p2_min = -Nmax + 1;  
	p2_max = Nmax  ;*/

    snprintf(drfilename1,128*sizeof(char), "ResFiles_N%ld/data_files/dataPlaneRing_%ld_%d_%d_R%d",NN, NN,p1,r1,NRi);
    drfile1= fopen(drfilename1, "rb"); 
    fread(buffer1, sizeof(double), 2 * NN*NN, drfile1);
    fclose(drfile1); 

    //for(p2=-Nmax+1;p2<Nmax;++p2) 
    for(p2=p2_min;p2<=p2_max;++p2)  
	{
	p3=-p1-p2; 
	sign_2=p2<0 ? -1 : 1 ;
	sign_3=p3<0 ? -1 : 1 ;	
	if(sign_3*p3<Nmax){

	   	double K[3][2]={{sqrt(k_perp_avg[r1]*k_perp_avg[r1] + p1*p1),p1},{sqrt(k_perp_avg[r2]*k_perp_avg[r2] + p2*p2),p2},{sqrt(k_perp_avg[r3]*k_perp_avg[r3] + p3*p3),p3}};
	   	
	   	sort(K);

// ###########  Parameterization for shape, size and orientation wrt LoS of triangle ################	
	if ( K[0][0]>(LL_k1) && K[0][0]<UL_k1){     	
	     	K1 = K[0][0];
	     	T  = K[1][0]/K[0][0];
	     	MU = 0.5*(T+1./T-(K[2][0]/K[0][0])*(K[2][0]/K[1][0]));
	     	pz = K[0][1]/K[0][0];
		px = ((K[1][1]/K[1][0]) + MU*pz)/(sqrt(1.-(MU*MU)));
		Num=rdata [rp * 4 + 3] ;

// ###########  Analytical redshift space bispectrum ####################
                mu_1 = K[0][1]/K[0][0] ;
                mu_2 = K[1][1]/K[1][0] ;
                mu_3 = K[2][1]/K[2][0] ;  	
	     	Biana=2.*ff*( input_ps(Cunits*K[0][0])*input_ps(Cunits*K[1][0]) + input_ps(Cunits*K[1][0])*input_ps(Cunits*K[2][0]) + input_ps(Cunits*K[2][0])*input_ps(Cunits*K[0][0]) );
	     	Biana = Biana*(1+beta*mu_1*mu_1)*(1+beta*mu_2*mu_2)*(1+beta*mu_3*mu_3);

if (MU==1)
{
printf("%d\t%d\t%d\t%d\t%d\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",r1,r2,r3,p1,p2,p3,K[0][0],K[1][0],K[2][0],MU,T,px,pz);
}
	     	
// ###########  Calculate bispectrum from data #########################	     	  
    		//FBI(p1,r1,p2,r2,p3,r3,sign_2,sign_3, &Bispec_real, &Bispec_imag);  
    		//FBI_fast(p1,r1,p2,r2,p3,r3,sign_2,sign_3, &Bispec_real, &Bispec_imag);  
    		//FBI_fast_new(p1,r1,p2,r2,p3,r3,sign_2,sign_3,buffer1, &Bispec_real, &Bispec_imag);

    		snprintf(drfilename2,128*sizeof(char), "ResFiles_N%ld/data_files/dataPlaneRing_%ld_%d_%d_R%d",NN, NN,sign_2*p2,r2,NRi);
		drfile2= fopen(drfilename2, "rb");
		snprintf(drfilename3,128*sizeof(char), "ResFiles_N%ld/data_files/dataPlaneRing_%ld_%d_%d_R%d",NN, NN,sign_3*p3,r3,NRi);
		drfile3= fopen(drfilename3, "rb");

		fread(buffer2, sizeof(double), 2 * NN*NN, drfile2);
		fread(buffer3, sizeof(double), 2 * NN*NN, drfile3);
		
		fclose(drfile2);
		fclose(drfile3);
		
		Bispec_real=0.; 
		Bispec_imag=0.;

		for(ii=0;ii<NN;++ii){
		  for(jj=0;jj<NN;++jj){
		      index=ii*(NN)+ jj;
		      Bispec_real+=buffer1[2 * index]*buffer2[2 * index]*buffer3[2 * index]-buffer1[2 * index]*sign_2*buffer2[2 * index + 1]*sign_3*buffer3[2 * index + 1]-buffer2[2 * index]*buffer1[2 * index + 1]*sign_3*buffer3[2 * index + 1]-buffer3[2 * index]*buffer1[2 * index + 1]*sign_2*buffer2[2 * index + 1];
				
		      Bispec_imag+=-buffer1[2 * index + 1]*sign_2*buffer2[2 * index + 1]*sign_3*buffer3[2 * index + 1]+buffer1[2 * index]*buffer2[2 * index]*sign_3*buffer3[2 * index + 1]+buffer2[2 * index]*buffer1[2 * index + 1]*buffer3[2 * index]+buffer3[2 * index]*buffer1[2 * index]*sign_2*buffer2[2 * index + 1];
		      

					}}

 //it give BS corresponding to particular triangle(shape,size,orintation wrt LoS) from data. At this point we have k1, mu, t , px, pz, Num, BS_real. Now we need to calculate multipoles of this BS_real and store them after binning

//printf("%d\t%d\t%d\t%d\t%d\t%d\t%lf\t%lf\t%lf\n",r1,r2,r3,p1,p2,p3,Bispec_real,Biana,Bispec_imag);
/*if (K1>Nmax-5 && K1<Nmax){  		 
printf("%d\t%d\t%d\t%d\t%d\t%d\t%lf\t%lf\t%lf\n",r1,r2,r3,p1,p2,p3,Bispec_real,Biana,Bispec_imag);
fp1=fopen(file1,"a");  
fprintf(fp1,"%d\t%d\t%d\t%d\t%d\t%d\t%lf\t%lf\t%lf\n",r1,r2,r3,p1,p2,p3,Bispec_real,Biana,Bispec_imag);
fclose(fp1);	
}*/
// fwrite(&K1, sizeof(double), 1, fk1);  // 35GB file for N=384

// ############# Binning  #######################################################

		K1_bin=floor((K1-(LL_k1))/(1.*BW_k1));
		T_bin  = (T==1) ? (n_mtbin-1) : floor ((T-.5)/(.5/(1.*n_mtbin)));
		MU_bin = (MU==1) ? (n_mtbin-1) : floor ((MU-.5)/(.5/(1.*n_mtbin)));
		
		//Bispec_real /= (vol*Num);
		//Bispec_imag /= (vol*Num);
		//Num /= (NN*NN);
		//printf("%lf\t\t\t%lf\t\t\t%f\n",Bispec_real,Biana,Num);
			
		cal_coeff(px, pz, bs, norm);
		

		if (Bispec_real==0.)// because we put \Delta(k=0)=0 by hand to avoid k=0 mode
		    {Biana=0.;
		     Num = 0.;}
		for(ii=0;ii<npoles;++ii)
		{
		bin_BS_real[K1_bin][T_bin][MU_bin][ii] += Bispec_real * bs[ii];
		bin_Biana[K1_bin][T_bin][MU_bin][ii] += Biana * Num * bs[ii];
		bin_normalization[K1_bin][T_bin][MU_bin][ii] += Num * norm[ii] ;
		}		} // if k1<Nmax loop

}}}
/*
if (rp%(rcount/Nmax)==0){
  printf("loop %ld out of %ld in  time %f\n",rp/(rcount/Nmax),Nmax,  omp_get_wtime()-time); 
  fp2=fopen(file2,"a");
  fprintf(fp2,"loop %ld out of %ld in  time %f\n",rp/(rcount/Nmax),Nmax,  omp_get_wtime()-time); 
  fclose(fp2);	
  time=omp_get_wtime();}   */ 

  printf("loop %d out of %ld in  time %f\n",p1,Nmax,  omp_get_wtime()-time); 
  fp2=fopen(file2,"a");
  fprintf(fp2,"loop %d out of %ld in  time %f\n",p1,Nmax,  omp_get_wtime()-time); 
  fclose(fp2);	
  time=omp_get_wtime();
}
printf("Total time %f\n", omp_get_wtime()-start_time); 		
fp2=fopen(file2,"a");
fprintf(fp2,"Total time %f\n", omp_get_wtime()-start_time); 
fclose(fp2);		 
  
  
  
  

for(ii=0;ii<npoles;++ii)
{  snprintf(file13, 128*sizeof(char),"ResFiles_N%ld/BS_test_ii%d_N%ld_R%d.txt",NN,ii,NN,NRi);
   fp13=fopen(file13,"w");
   
   for (i=0;i<n_k1bin;i++){
    for (j=0;j<n_mtbin;j++){
     for (l=0;l<n_mtbin;l++){
	bin_BS_real[i][j][l][ii]=bin_BS_real[i][j][l][ii]/(vol*bin_normalization[i][j][l][ii]);
	bin_Biana[i][j][l][ii]=bin_Biana[i][j][l][ii]/(bin_normalization[i][j][l][ii]);
	if (bin_normalization[i][j][l][0]>.1){
	  fprintf(fp13,"%d\t%d\t%d\t%e\t%e\t%f\n",i,j,l,bin_BS_real[i][j][l][ii],bin_Biana[i][j][l][ii],(bin_normalization[i][j][l][ii])/(NN*NN*1.));}
			}}} 
		
  fclose(fp13);
}
 
 
 
 
  

  
/*
for(ell=0;ell<=2;ell=ell+2)
{//ell=0;emm=0;
  norm_fact = 1.;//sqrt((2*ell+1)/ (4*pi));
  for(emm=0;emm<=ell;++emm)
	{//ell=2; emm=2; norm_fact = sqrt((2*ell+1)/ (4*pi));
	  ii= ell==2 ? ell+emm-1 : ell+emm;
	snprintf(file13, 128*sizeof(char),"BS_for_l=%d_m=%d.txt",ell,emm);
	fp13=fopen(file13,"w");
   for (i=0;i<Nmax;i++){
	for (j=0;j<n_mtbin;j++){
		for (l=0;l<n_mtbin;l++){
		bin_BS_real[i][j][l][ii]=bin_BS_real[i][j][l][ii]/(bin_normalization[i][j][l][ii]);
		bin_Biana[i][j][l][ii]=bin_Biana[i][j][l][ii]/(bin_normalization[i][j][l][ii]);
		if (bin_normalization[i][j][l][0]>.1)
		{
			fprintf(fp13,"%d\t%d\t%d\t%e\t%e\t%f\n",i,j,l,norm_fact*bin_BS_real[i][j][l][ii],norm_fact*bin_Biana[i][j][l][ii],(bin_normalization[i][j][l][ii]));
		}}}} 
		
		fclose(fp13);
}}*/






}






/*
for(ell=0;ell<=2;ell=ell+2){
    for(emm=0;emm<=ell;++emm)
	{
	  ii= ell==2 ? ell+emm-1 : ell+emm;
	  Ylm=gsl_sf_legendre_sphPlm(ell, emm, pz)*Expo(emm,px,pz);
	  bin_BS_real[K1_bin][T_bin][MU_bin][ii] += Bispec_real * Ylm;
	  bin_BS_imag[K1_bin][T_bin][MU_bin][ii] += Bispec_imag * Ylm;
	  bin_normalization[K1_bin][T_bin][MU_bin][ii] += Num*Ylm*Ylm;
	}}
*/

































