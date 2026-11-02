#include<stdlib.h>
#include<stdio.h>
#include<fftw3.h>
#include<math.h>
#include<unistd.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<omp.h>


extern long NN, Nr, Nmax;
extern int NRi;
extern double ***ro, ***RR,**D1, **D2, **D3, LL;
extern fftw_plan p;


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



double  input_ps(double kk)
// return input brightness temperature PS in K^2
{
  double pu;  

  double Av=1.0, betav=-2.0; // amplitude and exponent of PS 
  pu=Av*pow(kk,betav);
  return(pu);
}




/*-------------------  Df ----------------------*/

void Df(int plane,int ring, int sign, double **DD, fftw_plan pd)   
  	{
  	int i, j, index_1, index_2;
  	fftw_complex *to_DD, *comp_ro;
	to_DD=(fftw_complex*)&(DD[0][0]);
	comp_ro=(fftw_complex*)&(ro[0][0][0]);
	for (i=0;i<NN;++i)
  		for (j=0;j<NN;++j)
  		{	
			index_1=i*NN*(NN/2+1) + j*(NN/2+1)  + sign*plane;
			index_2=i*NN+j;
			to_DD[index_2][0]=comp_ro[index_1][0]*RR[ring][i][j];
			to_DD[index_2][1]=comp_ro[index_1][1]*RR[ring][i][j];	
		}
	if (plane==0) // k=0 mode exlusion
	    {
		to_DD[0][0]=0.;
		to_DD[0][1]=0.;
	    }

	fftw_execute(pd);	     
	}
	
	
	
	void sort(double K[3][2])
	{   int i,j;
		double swap1,swap2;
		for (i = 0; i < 2; i++) {  
       			for (j = 0; j < 2-i; j++){  
					if(K[j][0]<K[j+1][0]){
			            swap1=K[j][0];
			            swap2=K[j][1];
			            K[j][0]=K[j+1][0];
			            K[j][1]=K[j+1][1];
			            K[j+1][0]=swap1;
					    K[j+1][1]=swap2;
					} } }
	}
	
//	void sort(double *K)
//	{   int i,j;
//		double swap;
//		for (i = 0; i < 2; i++) {  
//       			for (j = 0; j < 2-i; j++){  
//					if(K[j]<K[j+1]){
//			            swap=K[j];
//			            K[j]=K[j+1];
//			            K[j+1]=swap;} } }
//	}
	
	
	
	
	
	
	double Expo(int m, double px, double pz)
	{
		double res, py;
		
		if(m==0)
		{
			res=1.;
		}
		
		if(m==1)
		{
			res=px/(sqrt(1-pz*pz));
		}
		
		if(m==2)
		{
			res=(pz*pz+2*px*px-1)/(1-pz*pz);
		}
		
		if(m==3)
		{
			res=(4*px*px*px + 3*px*pz*pz -3*px)/(pow(1-pz*pz,1.5));
		}
		
		if(m==4)
		{
			res=(8*pow(px,4)+pow(pz,4)+8*px*px*pz*pz-8*px*px-2*pz*pz+1)/(pow(1-pz*pz,2));
		}
		
		if(m==5)
		{
			res=(16*pow(px,5) + 20*pow(px,3)*pz*pz - 20*pow(px,3) - 10*px*pz*pz + 5*px + 5*px*pow(pz,4))/(pow(1-pz*pz,2.5));
		}
		
		if(m==6)
		{
			py=sqrt(1-px*px - pz*pz);
			res= (pow(px,6) - pow(py,6) - 15*pow(px,4)*py*py + 15*pow(py,4)*px*px)/(pow(1-pz*pz,3));
			//res=(36*pow(px,6)+ 48*pow(px,4)*pz*pz - 48*pow(px,4) - 36*pow(px,2)*pow(pz,2) + 18*px*px + 18*pow(pz,4)*px*px +3*pz*pz - 1 + pow(pz,6) - 3*pow(pz,4))/(pow(1-pz*pz,3));
		}
		
		return res;
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
//	
//	
//  	double **Df2(int plane,int ring, int sign)   
//  	{
//  	int i, j, index_1, index_2;
//  	
//  	double **BB;
//  	BB=allocate_2d_double(NN,NN*2); 
//	fftw_complex *to_BB, *comp_ro;
//	to_BB=(fftw_complex*)&(BB[0][0]);
//	comp_ro=(fftw_complex*)&(ro[0][0][0]);
//	fftw_plan q;
//  	q= fftw_plan_dft_2d(NN, NN, (fftw_complex*)&(BB[0][0]), (fftw_complex*)&(BB[0][0]),FFTW_BACKWARD, FFTW_ESTIMATE);
//  	
//  	
//	for (i=0;i<NN;++i)
//  		for (j=0;j<NN;++j)
//  		{	
//			index_1=i*NN*(NN/2+1) + j*(NN/2+1)  + sign*plane;
//			index_2=i*NN+j;
//			
//			
//			
////			to_AA[index_2][0]=comp_ro[index_1][0]*RR[ring][i][j];
////			to_AA[index_2][1]=comp_ro[index_1][1]*RR[ring][i][j];
//			to_BB[index_2][0]=comp_ro[index_1][0]*RR[ring][i][j];
//			to_BB[index_2][1]=comp_ro[index_1][1]*RR[ring][i][j];
//		}
//	//	 
//	//fftw_execute(q);   
//	    return BB;
//	}	
	
	
	
void fill_RR(double ***RR, double *k_perp_Low, double *k_perp_High) {
	int r,i,j,ia,ja;
	double kk;
    for ( r = 0; r < Nr; ++r) {
        for ( i = 0; i < NN; ++i) {
            for ( j = 0; j < NN; ++j) {
                 ia = (i > NN / 2) ? (i - NN) : i;
                 ja = (j > NN / 2) ? (j - NN) : j;
                 kk = sqrt(1. * (ia * ia + ja * ja));
                if (kk >= k_perp_Low[r] && kk < k_perp_High[r]) {
                    RR[r][i][j] = 1.0;
                } else {
                    RR[r][i][j] = 0.0;
                }
            }
        }
    }
}	


void read_field(double ***ro)
{	int i,j,k;	
	FILE *fp1;
	char file1[128];
  	snprintf(file1,128*sizeof(char),"ResFiles_N%ld/data_files/field_N%ld_R%d.txt",NN,NN,NRi);
   	fp1=fopen(file1,"r");
    	for(i=0;i<NN;++i)
      	  for(j=0;j<NN;++j)
            for(k=0;k<NN+2;++k){
		fread(&ro[i][j][k],sizeof(double),1,fp1); }
        
    	fclose(fp1);
}


void fill_sd(double sd[Nmax][Nr][2*NN*NN], int read_or_cal)
{
int r1,p1,ii,jj,index;
if (read_or_cal==0){
char drfilename1[128];
FILE *drfile1;
double buffer1[2*NN*NN];
printf("\n buffer1 defined safely\n");
printf("\n sd defined safely\n");
for(p1=0;p1<Nmax;++p1){
    for(r1=0;r1<Nr;++r1){
    snprintf(drfilename1,128*sizeof(char), "ResFiles_N%ld/data_files/dataPlaneRing_%ld_%d_%d_R%d",NN, NN,p1,r1,NRi);
    drfile1= fopen(drfilename1, "rb"); 
    fread(buffer1, sizeof(double), 2 * NN*NN, drfile1);
    fclose(drfile1);
for(ii=0;ii<NN;++ii){
  for(jj=0;jj<NN;++jj){
	index=ii*(NN)+ jj;
	sd[p1][r1][2 * index]=buffer1[2 * index];
        sd[p1][r1][2 * index + 1]=buffer1[2 * index + 1];
	}}
}}
}

else if (read_or_cal==1){
double **D;
fftw_plan pd;
D=allocate_2d_double(NN,NN*2); 
pd= fftw_plan_dft_2d(NN, NN, (fftw_complex*)&(D[0][0]), (fftw_complex*)&(D[0][0]),FFTW_BACKWARD, FFTW_ESTIMATE);
fftw_complex *to_D; 
to_D=(fftw_complex*)&(D[0][0]);


for(p1=0;p1<Nmax;++p1){
    for(r1=0;r1<Nr;++r1){
Df(p1,r1,1,D,pd); 
for(ii=0;ii<NN;++ii){
  for(jj=0;jj<NN;++jj){
	index=ii*(NN)+ jj;
	sd[p1][r1][2 * index] = to_D[index][0];
        sd[p1][r1][2 * index + 1] = to_D[index][1];
	}}
}}
fftw_destroy_plan(pd);
free(D);}
else{printf("Give either 0 to read or 1 to calculate\n");}
}	
	
	
// Writing all closed combinations of rings and corresponding triangles (multiplied by N^2).

void calCR(double  *k_perp_Low, double *k_perp_High, double *k_perp_avg)
{
int r,i,j,ia,ja,ii,jj,r1,r2,r3,index;
double  tt, mu,kk;
char filename[128];
snprintf(filename,128*sizeof(char), "ResFiles_N%ld/data_files/closed_rings_%ld", NN,NN);
FILE *file_rings = fopen(filename, "wb");
//FILE *file_rings  = fopen("closed_rings.dat", "w");
double data[4],Num;
int rcount;
rcount=0;

double ***II, **AA;
AA=allocate_2d_double(NN,NN*2); 
fftw_complex *to_AA;
to_AA=(fftw_complex*)&(AA[0][0]);
II= allocate_3d_double(Nr,NN,NN);



fftw_plan p;
p= fftw_plan_dft_2d(NN, NN, (fftw_complex*)&(AA[0][0]), (fftw_complex*)&(AA[0][0]),FFTW_BACKWARD, FFTW_ESTIMATE);	
	

	for (r=0;r<Nr;++r)
	  for(i=0;i<NN;++i)
	    for(j=0;j<NN;++j)
	    {	
	      	ia=(i>NN/2) ? (i-NN) : i ;
		ja=(j>NN/2) ? (j-NN) : j ; 
		kk=sqrt(1.*(ia*ia+ja*ja));
		if(kk>=k_perp_Low[r] && kk<k_perp_High[r])
			  {II[r][i][j]=1.;} 
		else  {II[r][i][j]=0.;}
			
		}


	
  	for (r=0;r<Nr;++r)
	{
	  for (i=0;i<NN;++i){
  		 for (j=0;j<NN;++j)
  			{
  				index=i*(NN)+ j;
  				to_AA[index][0]=II[r][i][j]; to_AA[index][1]=0.;
			}}
	fftw_execute(p);
	  for (i=0;i<NN;++i){
  		 for (j=0;j<NN;++j)
  			{
  				index=i*(NN)+ j;
  				II[r][i][j]=to_AA[index][0];
//  				printf("%e\n",to_AA[index][1]);
			}}}	
	





for(r1=0;r1<Nr;++r1){
 for(r2=0;r2<=r1;++r2){
	tt=k_perp_avg[r2]/k_perp_avg[r1];
	if(tt>=0.5 && tt<=1.){
	 for(r3=0;r3<=r2;++r3){
		mu=0.5*(tt+1./tt-(k_perp_avg[r3]/k_perp_avg[r1])*(k_perp_avg[r3]/k_perp_avg[r2]));
		if(mu>=0.5/tt && mu<=1.){
		Num=0.;
		for(ii=0;ii<NN;++ii){
	       	 for(jj=0;jj<NN;++jj)
		  {Num+=II[r1][ii][jj]*II[r2][ii][jj]*II[r3][ii][jj];}}
		
		//printf("%d\t%d\t%d\t%lf\n",r1,r2,r3,Num);
		data[0]=r1; data[1]=r2; data[2]=r3; data[3]=Num;
		fwrite(data, sizeof(double), 4, file_rings);
		rcount++ ;
		//fprintf(file_rings, "%lf %lf %lf %lf\n", data[0], data[1], data[2], data[3]);
	

		   }}}}} 
fclose(file_rings);	
	
char filename2[128];
snprintf(filename2,128*sizeof(char), "ResFiles_N%ld/data_files/rcountFor_%ld", NN,NN);
FILE *file_count = fopen(filename2, "w");
fprintf(file_count, "%d", rcount);
fclose(file_count);

printf("\n Total combinations of rings on a plane forming closed triangle = %d\n", rcount);



fftw_destroy_plan(p);
free(II);
free(AA);
}
	



void read_rings(int *rcountPtr, double **rdataPtr)	
{	
int i,j;	
// reading total number of confifigurations of closed rings 
char filename3[128];
snprintf(filename3,128*sizeof(char), "ResFiles_N%ld/data_files/rcountFor_%ld",NN, NN);
FILE *file_count = fopen(filename3, "r");
fscanf(file_count, "%d", rcountPtr);
fclose(file_count);

// reading data of total combinations
char rfilename[128];
snprintf(rfilename,128*sizeof(char), "ResFiles_N%ld/data_files/closed_rings_%ld", NN,NN);
FILE *rfile = fopen(rfilename, "rb");
*rdataPtr = (double *)malloc(*rcountPtr * 4 * sizeof(double));


for(i=0;i<*rcountPtr;++i){
for ( j = 0; j < 4; ++j) {
//fread(rdata[i], sizeof(double), 4, rfile);// this is required when data is read direcly in main code instead of this function
fread(*rdataPtr + i * 4 + j, sizeof(double), 1, rfile);}}
fclose(rfile);
//printf("\n I have read the rings combination data!!\n"); 
}	


/* // copy this in main code to implement
// writing i .dat file to compare with python
char filename[128];
snprintf(filename,128*sizeof(char), "ResFiles_N%ld/closed_rings_%ld.dat",NN, NN);
FILE *file_rings = fopen(filename, "wb");
   for (int i = 0; i < rcount; i++) {
        for (int j = 0; j < 4; j++) {
            fprintf(file_rings,"%lf\t ", rdata[i][j]);
        }
        fprintf(file_rings,"\n");}
fclose(file_rings);
*/
	
	
	



void wrapid()  // writes rings and planes ift of data
{
// writing IFT of each rings on each plane field in binary files i.e D

int p1,r1,ii,jj,index;
double dummy,dummy2;
double **D;
fftw_plan pd;
D=allocate_2d_double(NN,NN*2); 
pd= fftw_plan_dft_2d(NN, NN, (fftw_complex*)&(D[0][0]), (fftw_complex*)&(D[0][0]),FFTW_BACKWARD, FFTW_ESTIMATE);
fftw_complex *to_D; 
to_D=(fftw_complex*)&(D[0][0]);
char dfilename[128];

for(p1=0;p1<Nmax;++p1){
    for(r1=0;r1<Nr;++r1){
Df(p1,r1,1,D,pd); // note we are saving for +ve planes only. data for -ve planes can be achieved from these itself by negating imaginary part
snprintf(dfilename,128*sizeof(char), "ResFiles_N%ld/data_files/dataPlaneRing_%ld_%d_%d_R%d", NN,NN,p1,r1,NRi);
FILE *dfile= fopen(dfilename, "wb");
for(ii=0;ii<NN;++ii){
  for(jj=0;jj<NN;++jj){
	index=ii*(NN)+ jj;
	dummy=to_D[index][0];
	dummy2=to_D[index][1];
	//printf("%lf ---",dummy);
	fwrite(&dummy, sizeof(double), 1, dfile);
	fwrite(&dummy2, sizeof(double), 1, dfile);
	}}
fclose(dfile);}}
fftw_destroy_plan(pd);
free(D);
}

////////  printing outputs of wrapid	
/*
char drfilename3[128];
FILE *drfile3;
snprintf(drfilename3,128*sizeof(char), "ResFiles_N%ld/dataPlaneRing_%ld_%d_%d",NN, NN,1,1);
drfile3= fopen(drfilename3, "rb");
double buffer1[2*NN*NN];
fread(buffer1, sizeof(double), 2 * NN*NN, drfile3);
for(ii=0;ii<NN;++ii){
 for(jj=0;jj<NN;++jj){
	index=ii*(NN)+ jj;
	to_D2[index][0] = buffer1[2 * index];
	to_D2[index][1] = buffer1[2 * index + 1];}}
	///////////////////////
for(ii=0;ii<NN;++ii){
  for(jj=0;jj<NN;++jj){
	index=ii*(NN)+ jj;
	
	fread(&dummy, sizeof(double), 1, drfile3);
	fread(&dummy2, sizeof(double), 1, drfile3);
	to_D2[index][0]=dummy ;
	to_D2[index][1]=dummy2 ;}}
for(ii=0;ii<NN;++ii){
  for(jj=0;jj<NN;++jj){
	index=ii*(NN)+ jj;
	printf("%lf---",to_D2[index][0]);printf("%lf---",to_D2[index][1]);}
printf("\n");}
*/	



void wrapid_fast()  // writes rings and planes ift of data
{
// writing IFT of each rings on each plane field in binary files i.e D

int p1,r1,ii,jj,index;
double **D;
double buffer[2*NN*NN];
fftw_plan pd;
D=allocate_2d_double(NN,NN*2); 
pd= fftw_plan_dft_2d(NN, NN, (fftw_complex*)&(D[0][0]), (fftw_complex*)&(D[0][0]),FFTW_BACKWARD, FFTW_ESTIMATE);
fftw_complex *to_D; 
to_D=(fftw_complex*)&(D[0][0]);
char dfilename[128];

for(p1=0;p1<Nmax;++p1){
    for(r1=0;r1<Nr;++r1){
Df(p1,r1,1,D,pd); // note we are saving for +ve planes only. data for -ve planes can be achieved from these itself by negating imaginary part
snprintf(dfilename,128*sizeof(char), "ResFiles_N%ld/data_files/dataPlaneRing_%ld_%d_%d_R%d", NN,NN,p1,r1,NRi);
FILE *dfile= fopen(dfilename, "wb");
for(ii=0;ii<NN;++ii){
  for(jj=0;jj<NN;++jj){
	index=ii*(NN)+ jj;
	buffer[2 * index] = to_D[index][0];
        buffer[2 * index + 1] = to_D[index][1];
	}}
fwrite(buffer, sizeof(double), 2 * NN*NN, dfile);
fclose(dfile);}}
fftw_destroy_plan(pd);
free(D);
}




void FBI(int p1,int r1,int p2,int r2,int p3,int r3,int sign_2,int sign_3,double *Bispec_real, double *Bispec_imag)
	   	{
	   	double dummy_data[2];
    		int ii, jj, index;
    		
    		fftw_complex *to_D1, *to_D2, *to_D3; 
	
		to_D1=(fftw_complex*)&(D1[0][0]);
		to_D2=(fftw_complex*)&(D2[0][0]);
		to_D3=(fftw_complex*)&(D3[0][0]);
	   	char drfilename1[128];
		char drfilename2[128];
		char drfilename3[128];
		FILE *drfile1,*drfile2,*drfile3;
	   	snprintf(drfilename1,128*sizeof(char), "ResFiles_N%ld/data_files/dataPlaneRing_%ld_%d_%d_R%d",NN, NN,p1,r1,NRi);
		drfile1= fopen(drfilename1, "rb");
		snprintf(drfilename2,128*sizeof(char), "ResFiles_N%ld/data_files/dataPlaneRing_%ld_%d_%d_R%d",NN, NN,sign_2*p2,r2,NRi);
		drfile2= fopen(drfilename2, "rb");
		snprintf(drfilename3,128*sizeof(char), "ResFiles_N%ld/data_files/dataPlaneRing_%ld_%d_%d_R%d",NN, NN,sign_3*p3,r3,NRi);
		drfile3= fopen(drfilename3, "rb");
		
		for(ii=0;ii<NN;++ii){
		  for(jj=0;jj<NN;++jj){
			index=ii*(NN)+ jj;
			
			//fread(&dummy, sizeof(double), 1, drfile1); fread(&dummy2, sizeof(double), 1, drfile1);
			fread(&dummy_data, sizeof(double), 2, drfile1);
			to_D1[index][0]=dummy_data[0]; //dummy ;
			to_D1[index][1]=dummy_data[1]; //dummy2 ;
			
			//fread(&dummy, sizeof(double), 1, drfile2); fread(&dummy2, sizeof(double), 1, drfile2);
			fread(&dummy_data, sizeof(double), 2, drfile2);
			to_D2[index][0]=dummy_data[0]; //dummy ;
			to_D2[index][1]=dummy_data[1]; //dummy2 ;
			
			//fread(&dummy, sizeof(double), 1, drfile3); fread(&dummy2, sizeof(double), 1, drfile3);
			fread(&dummy_data, sizeof(double), 2, drfile3);
			to_D3[index][0]=dummy_data[0]; //dummy ;
			to_D3[index][1]=dummy_data[1]; //dummy2 ;
			}}
		fclose(drfile1);
		fclose(drfile2);
		fclose(drfile3);
		
		*Bispec_real=0.; 
		*Bispec_imag=0.;

		for(ii=0;ii<NN;++ii){
		  for(jj=0;jj<NN;++jj){
		      index=ii*(NN)+ jj;
		      *Bispec_real+=to_D1[index][0]*to_D2[index][0]*to_D3[index][0]-to_D1[index][0]*sign_2*to_D2[index][1]*sign_3*to_D3[index][1]-to_D2[index][0]*to_D1[index][1]*sign_3*to_D3[index][1]-to_D3[index][0]*to_D1[index][1]*sign_2*to_D2[index][1];
				
		      *Bispec_imag+=-to_D1[index][1]*sign_2*to_D2[index][1]*sign_3*to_D3[index][1]+to_D1[index][0]*to_D2[index][0]*sign_3*to_D3[index][1]+to_D2[index][0]*to_D1[index][1]*to_D3[index][0]+to_D3[index][0]*to_D1[index][0]*sign_2*to_D2[index][1];
					}}
					}
					
					














void FBI_fast(int p1,int r1,int p2,int r2,int p3,int r3,int sign_2,int sign_3,double *Bispec_real, double *Bispec_imag)
	   	{
	   	double buffer1[2*NN*NN];
	   	double buffer2[2*NN*NN];
	   	double buffer3[2*NN*NN];
    		int ii, jj, index;

	   	char drfilename1[128];
		char drfilename2[128];
		char drfilename3[128];
		FILE *drfile1,*drfile2,*drfile3;
	   	snprintf(drfilename1,128*sizeof(char), "ResFiles_N%ld/data_files/dataPlaneRing_%ld_%d_%d_R%d",NN, NN,p1,r1,NRi);
		drfile1= fopen(drfilename1, "rb");
		snprintf(drfilename2,128*sizeof(char), "ResFiles_N%ld/data_files/dataPlaneRing_%ld_%d_%d_R%d",NN, NN,sign_2*p2,r2,NRi);
		drfile2= fopen(drfilename2, "rb");
		snprintf(drfilename3,128*sizeof(char), "ResFiles_N%ld/data_files/dataPlaneRing_%ld_%d_%d_R%d",NN, NN,sign_3*p3,r3,NRi);
		drfile3= fopen(drfilename3, "rb");
		
		fread(buffer1, sizeof(double), 2 * NN*NN, drfile1);
		fread(buffer2, sizeof(double), 2 * NN*NN, drfile2);
		fread(buffer3, sizeof(double), 2 * NN*NN, drfile3);
		
		fclose(drfile1);
		fclose(drfile2);
		fclose(drfile3);
		
		*Bispec_real=0.; 
		*Bispec_imag=0.;

		for(ii=0;ii<NN;++ii){
		  for(jj=0;jj<NN;++jj){
		      index=ii*(NN)+ jj;
		      *Bispec_real+=buffer1[2 * index]*buffer2[2 * index]*buffer3[2 * index]-buffer1[2 * index]*sign_2*buffer2[2 * index + 1]*sign_3*buffer3[2 * index + 1]-buffer2[2 * index]*buffer1[2 * index + 1]*sign_3*buffer3[2 * index + 1]-buffer3[2 * index]*buffer1[2 * index + 1]*sign_2*buffer2[2 * index + 1];
				
		      *Bispec_imag+=-buffer1[2 * index + 1]*sign_2*buffer2[2 * index + 1]*sign_3*buffer3[2 * index + 1]+buffer1[2 * index]*buffer2[2 * index]*sign_3*buffer3[2 * index + 1]+buffer2[2 * index]*buffer1[2 * index + 1]*buffer3[2 * index]+buffer3[2 * index]*buffer1[2 * index]*sign_2*buffer2[2 * index + 1];
					}}
					
					
					
					}




















					
					
					
					
					
					
void cal_coeff(double px, double pz, double *bs, double *norm)
{		
bs[0]=1.;
bs[1]=2.*(3.*pz*pz-1.);
bs[2]=-sqrt(2./3.)* px*pz;
bs[3]=2.*sqrt(0.6666666666666666)*(-1. + 2.*pow(px,2) + pow(pz,2));
bs[4]=8.*(3. - 30.*pow(pz,2) + 35.*pow(pz,4));
bs[5]=(-4.*px*pz*(-3. + 7.*pow(pz,2)))/sqrt(5.);
bs[6]=4.* sqrt(0.4)*(-1 + 2.*pow(px,2) + pow(pz,2))*(-1. + 7.*pow(pz,2));
bs[7]=(-4.* px*pz*(-3. + 4.*pow(px,2) + 3.*pow(pz,2)))/sqrt(35.);
bs[8]=8.* sqrt(0.05714285714285714)*(8.*pow(px,4) + 8.*pow(px,2)*(-1. + pow(pz,2)) + pow(-1. + pow(pz,2),2));
bs[9]=16.*(-5. + 105.*pow(pz,2) - 315.*pow(pz,4) + 231.*pow(pz,6));
bs[10]=-8.*sqrt(0.09523809523809523)*px*pz*(5. - 30.*pow(pz,2) + 33.*pow(pz,4));
bs[11]=(32.*(-1. + 2.*pow(px,2) +pow(pz,2))*(1. - 18.*pow(pz,2) + 33.*pow(pz,4)))/sqrt(105.0);
bs[12]=(-16.*px*pz*(-3. + 4.*pow(px,2) + 3.*pow(pz,2))*(-3. + 11.*pow(pz,2)))/sqrt(105.0);
bs[13]=(16.*sqrt(0.2857142857142857)*(-1. + 11.*pow(pz,2))* (8.*pow(px,4) + 8.*pow(px,2)*(-1. + pow(pz,2)) + pow(-1. + pow(pz,2),2)))/3.;
bs[14]=(-16* pz*(16*pow(px,5) + 20*pow(px,3)*(-1 + pow(pz,2)) + 5*px*pow(-1 + pow(pz,2),2)))/(3.*sqrt(77.0));
bs[15]=(32.*(32.*pow(px,6) + 48.*pow(px,4)*(-1. + pow(pz,2)) + 18.*pow(px,2)*pow(-1. + pow(pz,2),2) + pow(-1. + pow(pz,2),3.)))/sqrt(231.);



norm[0]=1.;
norm[1]=((3.*pz*pz-1.))*((3.*pz*pz-1.));
norm[2]=pz*pz*(1.-pz*pz);
norm[3]=pow(-1. + pow(pz,2),2);
norm[4]=pow(3. - 30.*pow(pz,2) + 35.*pow(pz,4),2);
norm[5]=(1. - pow(pz,2))*pow(-3.*pz + 7.*pow(pz,3),2);
norm[6]=pow(-1. + pow(pz,2),2)*pow(-1. + 7.*pow(pz,2),2);
norm[7]=pow(pz,2)*pow(1. - pow(pz,2),3);
norm[8]=pow(-1. + pow(pz,2),4);
norm[9]=pow(-5. + 105.*pow(pz,2) - 315.*pow(pz,4) + 231.*pow(pz,6),2);
norm[10]=(1. - pow(pz,2))*pow(5.*pz - 30.*pow(pz,3) + 33.*pow(pz,5),2);
norm[11]=pow(-1. + pow(pz,2),2)*pow(1. - 18.*pow(pz,2) + 33.*pow(pz,4),2);
norm[12]=pow(1. - pow(pz,2),3)*pow(-3.*pz + 11.*pow(pz,3),2);
norm[13]=pow(-1 + pow(pz,2),4)*pow(-1 + 11*pow(pz,2),2);
norm[14]=pow(pz,2)*pow(1 - pow(pz,2),5);
norm[15]=pow(-1. + pow(pz,2),6);
}				
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					
					

void FBI_fast_old(int p1,int r1,int p2,int r2,int p3,int r3,int sign_2,int sign_3,double *Bispec_real, double *Bispec_imag)
	   	{
	   	double buffer1[2*NN*NN];
	   	double buffer2[2*NN*NN];
	   	double buffer3[2*NN*NN];
    		int ii, jj, index;
    		
    		fftw_complex *to_D1, *to_D2, *to_D3; 
	
		to_D1=(fftw_complex*)&(D1[0][0]);
		to_D2=(fftw_complex*)&(D2[0][0]);
		to_D3=(fftw_complex*)&(D3[0][0]);
	   	char drfilename1[128];
		char drfilename2[128];
		char drfilename3[128];
		FILE *drfile1,*drfile2,*drfile3;
	   	snprintf(drfilename1,128*sizeof(char), "ResFiles_N%ld/dataPlaneRing_%ld_%d_%d",NN, NN,p1,r1);
		drfile1= fopen(drfilename1, "rb");
		snprintf(drfilename2,128*sizeof(char), "ResFiles_N%ld/dataPlaneRing_%ld_%d_%d",NN, NN,sign_2*p2,r2);
		drfile2= fopen(drfilename2, "rb");
		snprintf(drfilename3,128*sizeof(char), "ResFiles_N%ld/dataPlaneRing_%ld_%d_%d",NN, NN,sign_3*p3,r3);
		drfile3= fopen(drfilename3, "rb");
		
		
		fread(buffer1, sizeof(double), 2 * NN*NN, drfile1);
		fread(buffer2, sizeof(double), 2 * NN*NN, drfile2);
		fread(buffer3, sizeof(double), 2 * NN*NN, drfile3);
		
		for(ii=0;ii<NN;++ii){
		  for(jj=0;jj<NN;++jj){
			index=ii*(NN)+ jj;
	
			to_D1[index][0] = buffer1[2 * index];
			to_D1[index][1] = buffer1[2 * index + 1];
			
			to_D2[index][0] = buffer2[2 * index];
			to_D2[index][1] = buffer2[2 * index + 1];
			
			to_D3[index][0] = buffer3[2 * index];
			to_D3[index][1] = buffer3[2 * index + 1];
			}}
		fclose(drfile1);
		fclose(drfile2);
		fclose(drfile3);
		
		*Bispec_real=0.; 
		*Bispec_imag=0.;

		for(ii=0;ii<NN;++ii){
		  for(jj=0;jj<NN;++jj){
		      index=ii*(NN)+ jj;
		      *Bispec_real+=to_D1[index][0]*to_D2[index][0]*to_D3[index][0]-to_D1[index][0]*sign_2*to_D2[index][1]*sign_3*to_D3[index][1]-to_D2[index][0]*to_D1[index][1]*sign_3*to_D3[index][1]-to_D3[index][0]*to_D1[index][1]*sign_2*to_D2[index][1];
				
		      *Bispec_imag+=-to_D1[index][1]*sign_2*to_D2[index][1]*sign_3*to_D3[index][1]+to_D1[index][0]*to_D2[index][0]*sign_3*to_D3[index][1]+to_D2[index][0]*to_D1[index][1]*to_D3[index][0]+to_D3[index][0]*to_D1[index][0]*sign_2*to_D2[index][1];
					}}
					}
				
					
					
					
					
					









void calpow_mom(int Nbin,double* power, double* kmode, double* power_P2,double* power_P4, double* power_P6,double *no)
{ 
  
  /******************** TO FIND POWER SPECTRUM **************************/
  
  long i, j, k, a, b, c, d,N1,N2,N3;
  long index, index1, index2;
  
  double *no2, *no4, *no6;
  fftw_complex *comp_ro;
  
  float fac1, fac2, fac3, m, mu, P2, P4,P6, scale;
  double norml;

  double pi=M_PI; 
  
  N1=NN;N2=NN;N3=NN;
  norml=1./(1.*N1*N2*N3);
  
  /*************** TAKING FOURIER TRANSFORM OF RO. **************/
  
  //fftwf_plan p_ro_dum=fftwf_plan_dft_r2c_3d(N1, N2, N3, &(ro_dum[0][0][0]), (fftwf_complex*)&(ro_dum[0][0][0]),FFTW_ESTIMATE);
  
  //fftwf_plan q_ro_dum=fftwf_plan_dft_c2r_3d(N1, N2, N3, (fftwf_complex*)&(ro_dum[0][0][0]), &(ro_dum[0][0][0]),FFTW_ESTIMATE);
  
  
  //fftwf_execute(p_ro_dum);
  
  comp_ro=(fftw_complex *)&(ro[0][0][0]);
  
  /*********** TO FIND POWER SPECTRUM OF RO. **************/
  
  no2=calloc((size_t)Nbin,sizeof(double));
  no4=calloc((size_t)Nbin,sizeof(double));
  no6=calloc((size_t)Nbin,sizeof(double));
  
  fac1=1./(1.*N1*N1);
  fac2=1./(1.*N2*N2);
  fac3=1./(1.*N3*N3);
  
  /**************** BINNING POWER SPECTRA **********************/
  
  for(i=0;i<Nbin;++i)
    {
      kmode[i]=0.0;
      
      power[i]=0.0;
      power_P2[i]=0.0;
      power_P4[i]=0.0;
      power_P6[i]=0.0;
      
      no[i]=0.0;
      no2[i]=0.0;
      no4[i]=0.0;
      no6[i]=0.0;
    }
  
  //***********************************************************************
  
  scale=log10(0.5*N1)/Nbin;
  
  //***********************************************************************

  /*---------------------- BINNING POWER SPECTRA -------------------*/
  
  /*-------------------------- half lines ------------------------- */
  
  for(i=1;i<=N1/2;i++)
    for(j=0;j<=N2/2;j=j+N2/2)
      for(k=0;k<=N3/2;k=k+N3/2)
	{
	  index = i*N2*(N3/2+1) + j*(N3/2+1) + k;
	  
	  m = sqrt(fac1*i*i + fac2*j*j + fac3*k*k);	      
	  
	  mu = (1.0*k)/(N3*m);
	  
	  P2 = 0.5*(3.*mu*mu - 1.);
	  P4 = 0.125*(35.0*powf(mu,4.0)-30.0*mu*mu + 3);
	  P6 = 0.0625*(231.0*powf(mu,6.0)-315.0*powf(mu,4.0)+105.0*mu*mu-5 );
	  
	  d=(int)floorf(log10(m*N1)/scale);  // logarithmic bins
	  
	  if(d>=0 && d<Nbin)
	    {
	      power[d]+=(double)(((comp_ro[index][0]* comp_ro[index][0])+(comp_ro[index][1]* comp_ro[index][1]))); 
	      power_P2[d]+=(double)(((comp_ro[index][0]* comp_ro[index][0])+(comp_ro[index][1]* comp_ro[index][1]))*P2);
	      power_P4[d]+=(double)(((comp_ro[index][0]* comp_ro[index][0])+(comp_ro[index][1]* comp_ro[index][1]))*P4);
	      power_P6[d]+=(double)(((comp_ro[index][0]* comp_ro[index][0])+(comp_ro[index][1]* comp_ro[index][1]))*P6);
	      
	      kmode[d]+=(double)(m);
	      no[d]= no[d] + (double)(1.0);
	      no2[d]=no2[d] + (double)(P2*P2);
	      no4[d]=no4[d] + (double)(P4*P4);
	      no6[d]=no6[d] + (double)(P6*P6);
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
	      
	      mu = (1.0*c)/(N3*m);
	      
	      P2 = 0.5*(3.*mu*mu - 1.);
	      P4 = 0.125*(35.0*powf(mu,4.0)-30.0*mu*mu + 3);
              P6 = 0.0625*(231.0*powf(mu,6.0)-315.0*powf(mu,4.0)+105.0*mu*mu-5 );
	      
	      d=(int)floorf(log10(m*N1)/scale);//logarithmic bins
	      
	      if(d>=0 && d<Nbin)
		{
		  power[d]+=(double)(((comp_ro[index][0]* comp_ro[index][0])+(comp_ro[index][1]* comp_ro[index][1]))); 
		  power_P2[d]+=(double)(((comp_ro[index][0]* comp_ro[index][0])+(comp_ro[index][1]* comp_ro[index][1]))*P2);
		  power_P4[d]+=(double)(((comp_ro[index][0]* comp_ro[index][0])+(comp_ro[index][1]* comp_ro[index][1]))*P4);
                  power_P6[d]+=(double)(((comp_ro[index][0]* comp_ro[index][0])+(comp_ro[index][1]* comp_ro[index][1]))*P6);
		  
		  kmode[d]+=(double)(m);
		  no[d]= no[d] + (double)(1.0);
		  no2[d]=no2[d] + (double)(P2*P2);
		  no4[d]=no4[d] + (double)(P4*P4);
                  no6[d]=no6[d] + (double)(P6*P6);
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
	      
	      mu = (1.0*c)/(N3*m);
	      
	      P2 = 0.5*(3.*mu*mu - 1.);
	      P4 = 0.125*(35.0*powf(mu,4.0)-30.0*mu*mu + 3);
              P6 = 0.0625*(231.0*powf(mu,6.0)-315.0*powf(mu,4.0)+105.0*mu*mu-5 );
	      
	      d=(int)floorf(log10(m*N1)/scale);//logarithmic bins
	      
	      if(d>=0 && d<Nbin)
		{
		  power[d]+=(double)(((comp_ro[index][0]* comp_ro[index][0])+(comp_ro[index][1]* comp_ro[index][1]))); 
		  power_P2[d]+=(double)(((comp_ro[index][0]* comp_ro[index][0])+(comp_ro[index][1]* comp_ro[index][1]))*P2);
		  power_P4[d]+=(double)(((comp_ro[index][0]* comp_ro[index][0])+(comp_ro[index][1]* comp_ro[index][1]))*P4);
                  power_P6[d]+=(double)(((comp_ro[index][0]* comp_ro[index][0])+(comp_ro[index][1]* comp_ro[index][1]))*P6);
		  
		  kmode[d]+=(double)(m);
		  no[d]= no[d] + (double)(1.0);
		  no2[d]=no2[d] + (double)(P2*P2);
		  no4[d]=no4[d] + (double)(P4*P4);
                  no6[d]=no6[d] + (double)(P6*P6);
		}
	    } /* end k for loop */
	  
	}
    }
  
  //***********************************************************************
  //***********************************************************************
  
  for(i=0;i<Nbin;i++)
    {
      if (no[i]>0.0) 
	{
	  power[i] = norml*power[i]/(1.0*no[i]*pow(LL,3.));	  
	  kmode[i]=2.*pi*kmode[i]/(no[i]*LL);
	}
      if(no2[i]>0.0)
	power_P2[i] = norml*power_P2[i]/(1.0*no2[i]*pow(LL,3.));
      
      if(no4[i]>0.0)
	power_P4[i] = norml*power_P4[i]/(1.0*no4[i]*pow(LL,3.));

      if(no6[i]>0.0)
	power_P6[i] = norml*power_P6[i]/(1.0*no6[i]*pow(LL,3.));
    }
  
  //***********************************************************************
  
//  for(i=0;i<N1;i++)
//    {
//      index1 = i*N2*(N3/2+1) ;	  
//      for(j=0;j<N2;j++)
//	{
//	  index2=index1 + j*(N3/2+1) ;
//	  for(k=0;k<(N3/2+1);k++)
//	    {
//	      index=index2 + k;
//	      comp_ro[index][0]=comp_ro[index][0]*norml;
//	      comp_ro[index][1]=comp_ro[index][1]*norml;
//	    }
//	}
//    }
  
  /*  now convert the array back to real space */
  
  //***********************************************************************
  
//  fftwf_execute(q_ro_dum);
//  
//  fftwf_destroy_plan(p_ro_dum);
//  fftwf_destroy_plan(q_ro_dum); 
  
} /* end function */











				
					
					
