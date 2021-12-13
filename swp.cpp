#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "ran3.h"
#include "vec.hpp"
using namespace std;

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)
#define PI 3.141592653589793238462643383

double s(vec&, double&, double&, double&, double&, double&);
double gamma(double, vec&);
vec rdt(double, vec&);
vec psi(vec&, double, vec*, int*, int, double, double, double, double, double);
bool rul(vec&, double, double, double);

int main()
{
    int N=100000;
    int Ne=20;
    int seed=-934121523;
    int tr,fps,nframes,ti[Ne];
    double xa,ya,edge,t,sigma,k,omega,phi,lambda,nlambda,d,delt,vs,ve,alpha,a,b;
    vec rp[N],rpdt[N],rc[Ne];
    ofstream file;

    /* initialize parameters */
    xa=3./2.;
    ya=9./8.;
    edge=1.2;
    tr=7;
    fps=30;
    nframes=fps*tr;
    
    sigma=0.2;
    nlambda=10;
    lambda=xa/nlambda;
    k=2*PI/lambda;
    
    a=1.25;
    b=1.5;
    alpha=3.0;
    omega=1;
    /*d=2*sqrt((edge*edge*xa*xa)+(edge*edge*ya*ya));*/
    d=(a+b)*edge*xa;
    delt=(d*k)/(alpha*omega*nframes);
    phi=PI/2;
    t=0;
    
    vs=omega/k;
    ve=alpha*vs;
    /*------------------------*/
    
    /* initialize emitter coordinates */
    for(int i=0;i<Ne;i++)
    {
        ti[i]=i*(nframes/Ne);
        rc[i](-a*edge*xa+ve*ti[i]*delt,0);
    }
    /* rc[0](); */
    
    /* initialize array of points */
    for(int i=0;i<N;i++)
	{
        rp[i](edge*xa*(2.*ran3(&seed)-1.),edge*ya*(2.*ran3(&seed)-1));
        rpdt[i]();
	}
    
    
    string pts="points/points_";
    for(int i=1;i<=nframes;i++)
    {
        string pu,snum;
        stringstream ssnum;
        ssnum<<i;
        snum=ssnum.str();
        ssnum.str("");
        pu=pts+snum;

        for(int j=0;j<N;j++)
        {
            rpdt[j]=psi(rp[j],Ne,rc,ti,i,delt,sigma,k,omega,phi);
        }
        
        file.open(pu.c_str());
        for(int j=0;j<N;j++)
        {
            file<<rpdt[j];
            if(j<N-1) file<<",";
            file<<endl;
        }
        file.close();
    }
	
    return 0;
}




double s(vec& r, double& t, double& sigma, double& k, double& omega, double& phi)
{
    return(2.0*PI*sigma/k)*cos((k*r.norm())-(omega*t)+phi);
}

double gamma(double s, vec& r)
{
    return((s/r.norm())+1.);
}

vec rdt(double gamma, vec& r)
{
    return(gamma*r);
}

vec psi(vec& rp, double N, vec* rc, int* ti, int ii, double delt, double sigma, double k, double omega, double phi)
{
    int Nc=0;
    double t;
    vec r,rsum;
    rsum();
    for(int i=0;i<N;i++)
    {
        r();
        r=rp-rc[i];
        t=(ii-ti[i])*delt;
        if (rul(r,t,k,omega)&&ii>=ti[i])
        {
            rsum+=(rdt(gamma(s(r,t,sigma,k,omega,phi),r),r)+rc[i]);
            Nc++;
        }
    }
    if(rsum.zero()) rsum=rp;
    else rsum-=((Nc-1)*rp);
    return(rsum);
}

bool rul(vec& r, double t, double k, double omega)
{
    double d,dd,delta;
    delta=0.5;
    d=(omega*t)/k;
    dd=(2*PI*delta)/k;
    return((r.norm()>=(d-dd))&&(r.norm()<=(d+dd)));
}





/*--------------------------------------------------------------------------*/
float ran3(int *idum)
//int *idum;
{
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj,mk;
	int i,ii,k;

	if (*idum < 0 || iff == 0) {
		iff=1;
		mj=MSEED-(*idum < 0 ? -*idum : *idum);
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1;i<=54;i++) {
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=ma[ii];
		}
		for (k=1;k<=4;k++)
			for (i=1;i<=55;i++) {
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
		inext=0;
		inextp=31;
		*idum=1;
	}
	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext]=mj;
	return mj*FAC;
}

#undef PI
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

