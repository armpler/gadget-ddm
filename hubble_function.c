#include <stdio.h>
#include <math.h>
#include <malloc.h>
#include <mpi.h>

#include "allvars.h"
#include  "proto.h"





double acch=1e-5;

double hubble_function(double a)
{
	double hubble;
	hubble=All.Hubble*sqrt((All.OmegaDm+All.OmegaBaryon)*pow(a,-3)+All.OmegaR*pow(a,-4)+All.OmegaLambda);
	return hubble;
}

#ifdef DDM
//In the dacaying model, unit of all parameters have h0, which is hubble paremeter of non-decaying model
double Hubble(double a,double y1,double y2)
{
	double hubble;
	hubble=All.OldHubble * Omegainit.h / All.OldHubbleParam * sqrt((y1+Omegainit.Omegab)*pow(a,-3)+y2*pow(a,-4)+Omegainit.Omegalambda);
	return hubble;
}
double f1(double t,double y1,double y2)           //equation 1
{
	double lifetime;
    lifetime=1/(All.decaylifetime*(3.15576e16))*All.UnitTime_in_s/All.OldHubbleParam;
	return lifetime*y1/((t+1)*Hubble(1/(1+t),y1,y2));
}

double f2(double t,double y1,double y2)           //equation 2
{
	double lifetime;
    lifetime=1/(All.decaylifetime*(3.15576e16))*All.UnitTime_in_s/All.OldHubbleParam;
	return -lifetime*y1/(pow(1 + t, 2)*Hubble(1/(1+t),y1,y2));
}
void rk(double t0,double t,double* y1,double *y2)              
{
	int n,j,i,k;
	int ifail=0;
	double dt;
	double *xh;
	double *x8;
	double *vh;
	double *v8;
	x8 == NULL;
	v8 == NULL;
	for(j=1;j<=30;j++)
	{
		n=100*(pow(2,j-1));
		n=n+1;
		xh=(double *) malloc(sizeof(double)*n);
		vh=(double *) malloc(sizeof(double)*n);
		dt=(t-t0)/(n-1);
		*xh=*y1;
		*vh=*y2;
		for(i=1;i<n;i++)
		{
			double k11=f1(t0+dt*(i-1),*xh,*vh);
			double k12=f2(t0+dt*(i-1),*xh,*vh);
			double k21=f1(t0+dt*(i-1)+dt/2,*xh+dt/2*k11,*vh+dt/2*k12);
			double k22=f2(t0+dt*(i-1)+dt/2,*xh+dt/2*k11,*vh+dt/2*k12);
			double k31=f1(t0+dt*(i-1)+dt/2,*xh+dt/2*k21,*vh+dt/2*k22);
			double k32=f2(t0+dt*(i-1)+dt/2,*xh+dt/2*k21,*vh+dt/2*k22);
			double k41=f1(t0+dt*(i-1)+dt,*xh+dt*k31,*vh+dt*k32);
			double k42=f2(t0+dt*(i-1)+dt,*xh+dt*k31,*vh+dt*k32);
			double a1=*xh+dt*(k11+2*k21+2*k31+k41)/6;
			double a2=*vh+dt*(k12+2*k22+2*k32+k42)/6;
			xh++,vh++;
			*xh=a1;
			*vh=a2;
			//printf("%.9lf %d\n",a1,j);
			//printf("%.9lf %d\n",a2,j);
		}
		xh=xh-n+1;
		vh=vh-n+1;
		//printf("%lf %d\n",*vh,j);
		ifail=0;
		if(j==1) ifail=1;
		else
		{
			for(k=1;k<=(1+(n-1)/2);k++)
			{
				if(ifail==0)
				{
					if((*xh>acch)&&(*x8>acch)&&(fabs(*x8/(*xh)-1.0)>acch)) ifail=1;
					if((*vh>acch)&&(*v8>acch)&&(fabs(*v8/(*vh)-1.0)>acch)) ifail=1;
					//printf("ifail=%d\n",ifail);
				}
				if(ifail==1)
				{
					x8=x8-k+1;
					v8=v8-k+1;
					//printf("%lf\n check",*v8);
					xh=xh-2*(k-1);
					vh=vh-2*(k-1);
					free(x8);
					free(v8);
					break;
				}
				x8=x8+1;
				v8=v8+1;
				xh=xh+2;
				vh=vh+2;
			}
		}
		if(ifail==0)
		{
			v8=v8-(n-1)/2-1;
			x8=x8-(n-1)/2-1;
			xh=xh-2;
			vh=vh-2;
			*y1=*xh;
			*y2=*vh;
			xh=xh-n+1;
			vh=vh-n+1;
			//printf("y1=%lf\n",*y1);
			//printf("y2=%lf\n",*y2);
			free(xh);
			free(x8);
			free(vh);
			free(v8);
			break;
		}
		x8=xh;
		v8=vh;
		//printf("%lf\n",*v8);
	}
}

void findh()
{
	double h0,h1,h2;
	h0=All.HubbleParam;
    All.OldHubbleParam = h0;
	int i,hstage;
	hstage=1;
	double omegadmfh,omegarfh,sum;
    double omegalambda0,omegab0,omegadm0,omegar0;
    if(ThisTask==0)
    {
    omegadm0=All.OmegaDm*h0*h0;
    omegalambda0=All.OmegaLambda*h0*h0;
    omegab0=All.OmegaBaryon*h0*h0;
    omegar0=All.OmegaR*h0*h0;
	sum = 0;
	for(i=1;fabs(sum-1)>(1e-10);i++)
	{
		if(hstage==1)
		{
			Omegainit.h=h0-0.1*(i-1);
			Omegainit.Omegalambda=omegalambda0/(Omegainit.h*Omegainit.h);
			Omegainit.Omegab=omegab0/(Omegainit.h*Omegainit.h);
			omegadmfh=omegadm0/(Omegainit.h*Omegainit.h);
			omegarfh=omegar0/(Omegainit.h*Omegainit.h);
			rk(1/All.TimeBegin-1,1/All.TimeMax-1,&omegadmfh,&omegarfh);
			sum=Omegainit.Omegalambda+Omegainit.Omegab+omegadmfh+omegarfh;
			if(sum>=1)
			{
				hstage=2;
				h1=Omegainit.h;
				h2=Omegainit.h+0.1;
			}
		}
		if(hstage==2)
		{
			Omegainit.h=(h1+h2)/2;
			Omegainit.Omegalambda=omegalambda0/(Omegainit.h*Omegainit.h);
			Omegainit.Omegab=omegab0/(Omegainit.h*Omegainit.h);
			omegadmfh=omegadm0/(Omegainit.h*Omegainit.h);
			omegarfh=omegar0/(Omegainit.h*Omegainit.h);
			rk(1/All.TimeBegin-1,1/All.TimeMax-1,&omegadmfh,&omegarfh);
			sum=Omegainit.Omegalambda+Omegainit.Omegab+omegadmfh+omegarfh;
			if(sum>=1) h1=Omegainit.h;
			if(sum<1)  h2=Omegainit.h;
			//printf("h=%lf\n",h);
		}
	}
	All.HubbleParam=Omegainit.h;
    All.OmegaBaryon=omegab0/(Omegainit.h*Omegainit.h);
    All.OmegaDm=omegadm0/(Omegainit.h*Omegainit.h);
    All.OmegaLambda=omegalambda0/(Omegainit.h*Omegainit.h);
    All.OmegaR=omegar0/(Omegainit.h*Omegainit.h);
    All.Hubble = All.OldHubble * Omegainit.h / All.OldHubbleParam;
    Omegainit.Omegab=omegab0/(Omegainit.h*Omegainit.h);
    Omegainit.Omegadm=omegadm0/(Omegainit.h*Omegainit.h);
    Omegainit.Omegalambda=omegalambda0/(Omegainit.h*Omegainit.h);
    Omegainit.Omegar=omegar0/(Omegainit.h*Omegainit.h);
    printf("h=%g\n",All.HubbleParam);
    printf("lifetime(internal units)=%g\n",1/(All.decaylifetime*(3.15576e16))*All.UnitTime_in_s/All.OldHubbleParam);
    }
    MPI_Bcast(&All, sizeof(All), MPI_BYTE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Omegainit,sizeof(Omegainit),MPI_BYTE,0,MPI_COMM_WORLD);

}

//function for gauss integral
double integral_function(double t0,double t)
{
    double dm1,r1;
    dm1=All.OmegaDm;
    r1=All.OmegaR;
    rk(t0,t,&dm1,&r1);
    return -1/((1+t)*Hubble(1/(1+t),dm1,r1));
}

double gauss_integral(double t0,double t1)  //5-point gauss integral
{
    double weight[5];
    double points[5];
    points[0]=0;
    weight[0]=128.0/225;
    points[1]=-1/3.0*sqrt(5-2*sqrt(10.0/7));
    weight[1]=(322+13*sqrt(70))/900;
    points[2]=1/3.0*sqrt(5-2*sqrt(10.0/7));
    weight[2]=(322+13*sqrt(70))/900;
    points[3]=-1/3.0*sqrt(5-2*sqrt(10.0/7));
    weight[3]=(322-13*sqrt(70))/900;
    points[4]=1/3.0*sqrt(5-2*sqrt(10.0/7));
    weight[4]=(322-13*sqrt(70))/900;
    int i,j,k,n;
    double y,dt;
    double ylast;
    double tl,tr;
    for(j=1;j<=30;j++)
    {
        n=100*pow(2,j-1);
        dt=(t1-t0)/n;
        y=0;
        for(i=1;i<=n;i++)
        {
            tl=t0+dt*(i-1);
            tr=tl+dt;
            for(k=0;k<=4;k++)
            {
                y=y+dt/2*weight[k]*integral_function(t0,dt/2*points[k]+(tl+tr)/2);
            }
        }        
	if(j==1) ylast=y;
        else
        {
            if(fabs(y/ylast-1)<1e-8) break;
            ylast=y;
        }
    }
    return y;
}
/* find omega of current time in decaying dark matter and reduce the mass of particles*/
void find_new_omega_and_change_mass()
{
    if(All.TimeStep > 1e-5)
    {
    if(ThisTask == 0) printf("finding new omega and changing mass\n");
    double newomegadm,newomegar;
    newomegadm=Omegainit.Omegadm;
    newomegar=Omegainit.Omegar;
    double reducefactor;
    double gamma;
    double Time_to_now;
    if(ThisTask==0)
    {
    gamma=1/(All.decaylifetime*(3.15576e16))*All.UnitTime_in_s   / All.OldHubbleParam;
    //printf("%.9f\n", All.Time);
    //printf("%.9f\n", All.TimeStep);
    Time_to_now=gauss_integral(1/(All.Time-All.TimeStep)-1,1/All.Time-1);
    printf("%.20lf\n", -gamma*Time_to_now);
    reducefactor=(1-All.OmegaDm/(All.OmegaDm+All.OmegaBaryon))+All.OmegaDm/(All.OmegaDm+All.OmegaBaryon)*exp(-gamma*Time_to_now);
    printf("reducefactor = %.20f\n", reducefactor);
    }
    MPI_Bcast(&reducefactor,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    int i;
    for(i=0;i<NumPart;i++)
    {
        P[i].Mass=P[i].Mass * reducefactor;
    }
    header.mass[1]=header.mass[1] * reducefactor;
    All.MassTable[1]=All.MassTable[1] * reducefactor;
    MPI_Barrier(MPI_COMM_WORLD);
    rk(1/All.TimeBegin-1,1/All.Time-1,&newomegadm,&newomegar);
    All.OmegaDm=newomegadm;
    All.OmegaR=newomegar;
    }
}
#endif
