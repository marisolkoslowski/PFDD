/*
 *
 *  Created by marisol 
 *  Copyright 2008 Purdue. All rights reserved.
 *
 */

#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#define N1 32 /*input number of grid points must be power of 2*/
#define N2 32 /*input number of grid points must be power of 2*/
#define N3 32 /*input number of grid points must be power of 2*/
#define NS 2  /*  number of slip systems for 12 for FCC */
#define ND 3
#define NT 100
#define NP 2  /*input number of stress increments*/
#define pi 3.141592654 

void fourn(float *, int nn[4], int, int);

int main(void)
{
int nsize, it, itp, isa, isb, psys, nf[4], i, j, k, na0, na, na1, nad1, nad2, nb;
double b, C11, C12, C44, S11, S12, S44, Asf[NS], Young, Poisson, Peierls, dslip, CD;
double L, d1, d2, d3, size, size3;
double setsigma, sigstep, setobs ,sigma[ND][ND], tau[NS];
double xn[NS][ND], xb[NS][ND], eps[NS][ND][ND];
double *fx, *fy, *fz;
double *BB , *xi, *xi_bc;
float *data, *data2;

void set3D2sys(double xn[NS][ND], double xb[NS][ND]);
void frec (double *, double *, double *, double, double,double);
void seteps(double eps[NS][ND][ND], double xn[NS][ND], double xb[NS][ND],double);
void Bmatrix( double *, double*, double*, double*, double eps[NS][ND][ND],  
			   double , double , double , double, double , double , double, double);
void initial(float *, double *, double *, double);	
void resolSS(double sigma[ND][ND], double tau[NS], double xn[NS][ND], double xb[NS][ND]);

    /* data arrays*/

  fx = malloc((N1)*(N2)*(N3)*sizeof(double));
  fy = malloc((N1)*(N2)*(N3)*sizeof(double));
  fz = malloc((N1)*(N2)*(N3)*sizeof(double));
  BB = malloc((NS)*(NS)*(N1)*(N2)*(N3)*sizeof(double));
  xi = malloc(2*(NS)*(N1)*(N2)*(N3)*sizeof(double));
  xi_bc = malloc(2*NS*N1*N2*N3*sizeof(double));

  data =  malloc((2*(NS)*(N1)*(N2)*(N3)+1)*sizeof(float));
  data2 =  malloc((2*(NS)*(N1)*(N2)*(N3)+1)*sizeof(float));

 /* data is offset by one respect to xi due to offset in fourn
     data goes to fft and is repalced by its fft
     xi is always in real space  and is updated every step*/

    /* Material constants isotropic Ni */
	
       b = 3.52E-10/sqrt(2.0); /* input Grid size = Burgers vector*/
	
    
	Young = 20E10; /* input Young modulus*/
	Poisson = 0.3; /*input Poisson ratio*/
	
		
	C44 = Young/2.0/(1.0+Poisson);
	C12 = 2.0 * C44 * Poisson/(1.0-2.0* Poisson);
	C11 = C12 + 2.0* C44;  
	S11 = (C12+C44)/C44/(3 * C12 + 2 * C44);
	S12 = -C12/2.0/C44/(3 * C12 + 2 * C44);
	S44 = 1.0/C44;

        dslip = 1.0 ;  /*in units of b*/
        CD = 1.0; /* kin coefficient for G-L equation not important for steady state*/ 
	Peierls = 0.6; /*input  in Joules Peierls energy */
          /* sinusoidal stacking fault energy coefficient*/
	Asf[0] = Peierls/(C44*dslip*b); 
	Asf[1] = Peierls/(C44*dslip*b);


    /* Setting grid parameters  */ 

	setobs = 0.30; //set random obstacle density

	L = (double)(N3-1);
	
	nf[0]=N1;
	nf[1]=N1;
	nf[2]=N2;
	nf[3]=N3;
	nsize = N1*N2*N3;
	
 	d1 = 1.0; //in unis of b so Fequencies are normalized
	d2 = 1.0;
	d3 = 1.0;
	
	size = L*b;
	size3 = size*size*size;

    /* Setting loading step */ 
   

	setsigma = 0.0000;   //initial stress
	sigstep = 0.002; /* input stress increment normalized with respect to shear modulus*/

       /*Set the slip systems */

	set3D2sys(xn,xb); 

        seteps(eps, xn, xb, dslip);

  	/* Frequencies for FFT */

 	 frec(fx, fy, fz, d1, d2, d3);


	 /* Setting the B matrix*/

        Bmatrix(BB, fx, fy, fz, eps, d1,  d2 , d3,  C11,  C12,  C44, b, dslip);

	

       printf(" Generating initial data\n");

       initial( data, xi, xi_bc, setobs);
 
      /* time step */

     printf("start loading step \n");  
 
     for(itp=0;itp<NP;itp++)
    {

	  /* itp = loading step —— Applied stress and resolved shear stress*/
	  
	  sigma[0][0]=0.0;
	  sigma[0][1]=0.0;
	  sigma[0][2]=1.0*(setsigma + itp*sigstep); 
	  sigma[1][0]= 0.0;
	  sigma[1][1]=0.0;
	  sigma[1][2]=0.0;
	  sigma[2][0]=sigma[0][2];
	  sigma[2][1]=0.0;
	  sigma[2][2]=0.0;

	  resolSS(sigma,tau, xn, xb);	


  /* it is the time step in G-L equation*/ 

      for(it=0;it<NT;it++)   
      {
      printf("loading step  - > time step, %d %d\n", itp, it);

	for(isa=0;isa<NS;isa++)
	  {
	  psys = 2*(isa*N1*N2*N3);
	  fourn(&data[psys],nf,3,-1);  /* FFT*/
	  }
	

      for (i=0; i<2*N1*N2*N3*NS+1; i++)
	 {
	  data2[i] = 0;
	 }

      for(isa=0;isa<NS;isa++)
	{
	  for(isb=0;isb<NS;isb++)
	    { 
	      for(i=0;i<N1;i++)
		for(j=0;j<N2;j++)
		  for(k=0;k<N3;k++)
		    {
		      na0 = 2*(i*N2*N3 + j*N3 + k + isb*N1*N2*N3);
		      na = 2*(i*N2*N3 + j*N3 + k + isa*N1*N2*N3);
		      na1 = na0+1;
		      nad1 = na0+1;
		      nad2 = na0+2;
		      nb = k+(j)*N3+(i)*N2*N3+(isa)*N1*N2*N3+(isb)*N1*N2*N3*NS;		
		      data2[na+1] += data[nad1] * BB[nb];
		      data2[na+2] += data[nad2] * BB[nb];
		    }
	    }
	 }



	for(isa=0;isa<NS;isa++)
	{
        for(i=0;i<N1;i++)
	    for(j=0;j<N2;j++)
	      for(k=0;k<N3;k++)
		{
		  na0 = 2*(k+(j)*N1+(i)*N1*N2+(isa)*N1*N2*N3);		 
		  na = 2*(k+(j)*N1+(i)*N1*N2);
		  na1 = na0+1;
		  nad1 = na0+1;
		  nad2 = na0+2;
		  if(xi_bc[na0]==0){
		  xi[na0] = xi[na0]-CD*(data2[nad1]/nsize+Asf[isa]*sin(2.0*pi*xi[na0])-tau[isa]/dslip);  
		  xi[na1] = xi[na1]-CD*(data2[nad2]/nsize);
		  } 
		  data[nad1] = xi[na0];
		  data[nad2] = xi[na1];
			  
		 } /* end ijk*/
	  } /* end isa */
		
     }  /*end it*/
   }  /*end itp*/
  
   return 0; 	 
}





 /* functions*/ 


#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void fourn(float data[], int nn[], int ndim, int isign)
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
		/*printf("%u, %u, %u, %u\n", ip1, ip2, ip3, nprev);*/
		for (i2=1;i2<=ip2;i2+=ip1) {
			if (i2 < i2rev) {
			  /*printf("%u, %u, %u, %u\n", i1, i2, i2rev, ibit);*/
				for (i1=i2;i1<=i2+ip1-2;i1+=2) {
				  /*printf("i1 %u, %u\n", i1, i2);*/
					for (i3=i1;i3<=ip3;i3+=ip2) {
						i3rev=i2rev+i3-i2;
						SWAP(data[i3],data[i3rev]);
						SWAP(data[i3+1],data[i3rev+1]);
					}
				}
			}
			ibit=ip2 >> 1;
			/*printf("ibit %u\n", ibit);*/
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
		/*printf("nprev %u, %u\n", nprev, n);*/
	}
}
void set3D2sys(double xn[NS][ND], double xb[NS][ND])
{
  xn[0][0]=0.0;
  xn[0][1]=0.0;
  xn[0][2]=1.0;
  
  xb[0][0]=0.5;
  xb[0][1]=sqrt(3.0)/2.0;
  xb[0][2]=0.0;
  
  xn[1][0]=0.0;
  xn[1][1]=0.0;
  xn[1][2]=1.0;
	
  xb[1][0]=0.5;
  xb[1][1]=-sqrt(3.0)/2.0;
  xb[1][2]=0.0;
  
  printf("Burgers vector b1 = (%lf , %lf , %lf )\n", xb[0][0],xb[0][1], xb[0][2]);
  printf("Burgers vector b2 = (%lf , %lf , %lf )\n", xb[1][0],xb[1][1], xb[1][2]);
  printf("Slip plane     n  = [%lf,  %lf , %lf ]\n", xn[0][0],xn[0][1], xn[0][2]);
  return;
}



void frec( double *fx,double *fy, 
		  double *fz, double d1, double d2, double d3)
{
	int i,j,k,ksym, nf;


for(i=0;i<N1;i++)
{
	for(j=0;j<N2;j++)
	{
		for(k=0;k<N3;k++)
		{
			nf = k+(j)*N3+(i)*N3*N2;
			/* frecuency in x */
			if (i==0) {
				fx[nf]= 0.0;
			}
			if (i >= 1 && i < N1/2 ) {
				fx[nf]= (double)(i)/(double)(N1)/d1;
			}
			if (i >= N1/2) {
				fx[nf]= ((double)(i)-(double)(N1))/(double)(N1)/d1;
				
			}
			/* frecuency in y */
			if (j==0) {
				fy[nf]= 0.0;
			}
			if (j >= 1 && j < N2/2 ) {
				fy[nf]= (double)(j)/(double)(N2)/d2;
			}
			if (j >= N2/2) {
				fy[nf]= ((double)(j)-(double)(N2))/(double)(N2)/d2;
			}				
			/* frecuency in z */
			if (k==0) {
				fz[nf]= 0.0;
			}
			if (k >= 1 && k < N3/2 ) {
				fz[nf]= (double)(k)/(double)(N3)/d3;
			}	
			if (k >= N2/2) {
				fz[nf]= ((double)(k)-(double)(N2))/(double)(N2)/d2;
			}		
			
				/*printf("%d %d %d    %lf %lf %lf \n", i, j, k, fx[nf], fy[nf],fz[nf]); 	*/

			}
		}
	}
	 
	
	return;
}


void seteps (double eps[NS][ND][ND], double xn[NS][ND], double xb[NS][ND],double dslip)
{
	int i, j, k;
	
	
	for (i=0; i<ND; i++) {
		for (j=0; j<ND; j++) {
			for (k=0; k<NS;k++){
				eps[k][i][j]= xb[k][i]*xn[k][j]/ dslip;				
			}
		}
	}
	
return;
}

void Bmatrix (double *BB, double *fx, double *fy, double *fz, double eps[NS][ND][ND],double d1, double d2,double d3, double C11, double C12, double C44, double b, double dslip)
{
printf("setting B matrix \n");
#define 	DELTA(i, j)   ((i==j)?1:0)
	
  int i,j,k,l,m,n, u, v, k1,k2,k3,ka,kb, nv, nb, nfreq;
  int is, js, ks;
  double  fkr;
  double C[ND][ND][ND][ND];
  float A[ND][ND][ND][ND];
  float B[NS][NS][N1][N2][N3];
  float G[ND][ND];
  double fk[ND];
  double xnu, fk2, fk4, fka,fkb;
  
  xnu = C12/(2.0*(C44+C12));
  
  /* set Cijkl*/
  
  for (i=0; i<ND; i++) {
    for (j=0; j<ND; j++) {
      for (k=0; k<ND; k++) {
	for (m=0; m<ND; m++) {
	  C[i][j][k][m] = C44 * (DELTA(i,k)*DELTA(j,m)+DELTA(i,m)*DELTA(j,k))+C12*DELTA(i,j)*DELTA(k,m);
	  A[i][j][k][m] = 0.0;
	  /*printf("in BB %d %d %d %d %d\n",ND, i,j,k,m);*/
	}
      }
    }
  }
  
  /* set A, Green function and B matrix*/
	
  for(k1=0;k1<N1;k1++)
    {
      for(k2=0;k2<N2;k2++)
	{
	  for(k3=0;k3<N3;k3++)
	    {
	      nfreq = k3+(k2)*N3+(k1)*N3*N2;
	      fk[0] = fx[nfreq];
	      fk[1] = fy[nfreq];
	      fk[2] = fz[nfreq];
	      fk2 = fk[0]*fk[0]+fk[1]*fk[1]+fk[2]*fk[2];
	      fk4 = fk2*fk2;
	      if(fk2>0)
		{
		  for (m=0; m<ND; m++) {
		    for (n=0; n<ND; n++) {
		      for (u=0; u<ND; u++) {
			for (v=0; v<ND; v++) {
			  A[m][n][u][v] = 0.0;
			  
			  for	(i=0; i<ND; i++) {
			    for (j=0; j<ND; j++) {
			      for (k=0; k<ND; k++) {
				G[k][i] = (2.0 * DELTA(i,k)/fk2-1.0/(1.0-xnu)*fk[i]*fk[k]/fk4)/(2.0*C44);
				for	(l=0; l<ND; l++) {
				  A[m][n][u][v] = A[m][n][u][v] - C[k][l][u][v]*C[i][j][m][n]*G[k][i]*fk[j]*fk[l] ;	
				}	
			      }
			    }
			  }
			  A[m][n][u][v] = A[m][n][u][v]+C[m][n][u][v];
			  /*printf("A %d %d %d %d %f %f %f %f\n",m,n,u,v,fk[0],fk[1],fk[2],A[m][n][u][v]);*/
			}
		      }
		    }
		  }
		  
		} /*if fk2 */										
	      for(ka=0;ka<NS;ka++)
		{
		  for(kb=0;kb<NS;kb++)	    
		    {
		      B[ka][kb][k1][k2][k3] = 0.0;
		      for (m=0; m<ND; m++) {
			for (n=0; n<ND; n++) {
			  for (u=0; u<ND; u++) {
			    for (v=0; v<ND; v++) {
				
			      B[ka][kb][k1][k2][k3]= B[ka][kb][k1][k2][k3] + A[m][n][u][v]*eps[ka][m][n]*eps[kb][u][v];
						      
			    }	
			  }
			}
		      }	  
		      nb = nfreq +(ka)*N1*N2*N3+(kb)*N1*N2*N3*NS;
		      BB[nb] = B[ka][kb][k1][k2][k3]/C44;
		      /*printf("%lf %lf %lf %lf \n", fx[nfreq], fy[nfreq], fz[nfreq], BB[nb]);*/
		    } /*kb*/
		}/* ka*/ 
	    }/*k3*/
	}/*k2*/	
    }/*k1*/
  return;
}

void initial(float * data, double * xi, double * xi_bc, double setobs)
{
	int ir, na0, na, na1,nad1,nad2, is, i, j, k, in;
	float r, rmin, rho;

	
	for(is=0;is<NS;is++)
    {
		printf("initial data for system %d\n",is);
		for(i=0;i<N1;i++)
		{	
			for(j=0;j<N2;j++)
			{
				for(k=0;k<N3;k++)
				{		
					na0 = 2*(k+(j)*N1+(i)*N1*N2+(is)*N1*N2*N3);
					na = 2*(k+(j)*N1+(i)*N1*N2);
					na1 = na0+1;
					nad1 = na0+1;
					nad2 = na0+2;
					
					xi[na0] = 0.0;      
					xi[na1] = 0.0;
					xi_bc[na0]=1.0; /* xi_bc=0 evolves  xi_bc=1 is fixed*/
						
						
					data[nad1] = xi[na0];
					data[nad2] = xi[na1];
					
				}/*k*/
			}/*j*/
		}/*i*/
			
    } /*is*/
	
	return;
}

void resolSS(double sigma[ND][ND], double tau[NS], double xn[NS][ND], double xb[NS][ND])
{
  int is, i, j, k;
  for (is=0;is<NS;is++){
    tau[is] = 0.0;
    for(i=0;i<ND;i++){
      for(j=0;j<ND;j++){
	tau[is] = tau[is]+sigma[i][j]*xn[is][j]*xb[is][i];	
	/*printf("%d %d %lf %lf %lf\n", i,j,sigma[i][j],xn[is][j], xb[is][i]);*/
      }	
    }
    /*printf("Resolved shear stresses\n”);*/
    printf("tau[%d] = %lf\n",is, tau[is]);
  }
  return;
}


