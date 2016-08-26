#ifndef RELAYLIB_H
#define RELAYLIB_H
#include<math.h>
#include<complex.h>
#define PI 3.141592653589793

void simetrica(double complex X[3], double complex S[3])
{
	double complex a=cos(2*PI/3)+sin(2*PI/3)*I;
	double complex A[3][3]={1,1,1,1,a*a,a,1,a,a*a};
	double complex invA[3][3]={1/3,1/3,1/3,1/3,a/3,a*a/3,1/3,a*a/3,a/3};
	S[0]=(X[0]+X[1]+X[2])/3;
	S[1]=(X[0]+a*X[1]+a*a*X[2])/3;
	S[2]=(X[0]+a*a*X[1]+a*X[2])/3;
}
double complex dft(double x[], double Fs, double Fk, int N)
{
	double ret_r=0;
	double ret_i=0;
	double complex RET;
	int k,i;
	k=(Fk*N/Fs);
	for(i=0; i<N;i++)
	{
		ret_r=2*x[i]*cos(2*PI*k/N*i)/N+ret_r;
		ret_i=-2*x[i]*sin(2*PI*k/N*i)/N+ret_i;
	}
	RET=ret_r+ret_i*I;
	return RET;
}
complex extractDC(complex X, double tau, int N, double dt)
{
	complex ret;
	double E1, F1, a1, fi1;
	E1=1-(1/exp(dt/tau))*cos(2*PI/N);
	F1=((1/exp(dt/tau)))*sin(2*PI/N);
	a1=sqrt(E1*E1+F1*F1);
	fi1=atan2(F1,E1);
	a1=1;
	fi1=0;
	ret=(1,1);
	return ret;
}
double RMS(double M[], int n)
{
	double ret;
	double sum=0;
	int i;
	for(i=0; i<n ; i++)
	{
	sum=M[i]*M[i]/n+sum;
	}
	ret=sqrt(sum);
	return ret;
}
double FaultDetector(double angA,double angB,double angC,double R0f,double R2f)
{
	double ret;
	ret=0;
	if((angA >= 0 && angA <= 60) && (angB >= 120 && angB <=180 )
	&& (angC >= 60 && angC <= 120 ) && (R0f >=0.1 && R0f <= 1.5 )
	&& (R2f>=0.22 && R2f <= 1.5))
	ret= 9;
	if((angA >= 60 && angA <= 120) && (angB >= 0 && angB <=60 )
	&& (angC >= 60 && angC <= 120 ) && (R0f >=0.1 && R0f <= 1.5 ) // ERROR DETECTADO
	&& (R2f>=0.22 && R2f <= 1.5))
	ret= 5;
	if((angA >= 120 && angA <= 180) && (angB >= 60 && angB <=120 )
	&& (angC >= 0 && angC <= 60 ) && (R0f >=0.1 && R0f <= 1.5 )
	&& (R2f>=0.22 && R2f <= 1.5))
	ret= 3;
	if((angA >= 0 && angA <= 60) && (angB >= 60 && angB <=120 )
	&& (angC >= 120 && angC <= 180 ) && (R0f >=0 && R0f <= 0.1 )
	&& (R2f>=0.22 && R2f <= 1.5))
	ret= 12;
	if((angA >= 120 && angA <= 180) && (angB >= 0 && angB <=60 )
	&& (angC >= 0 && angC <= 60 ) && (R0f >=0 && R0f <= 0.1 )
	&& (R2f>=0.71 && R2f <= 1.5))
	ret= 6;
	if((angA >= 60 && angA <= 120) && (angB >= 120 && angB <=180 )
	&& (angC >= 0 && angC <= 60 ) && (R0f >=0 && R0f <= 0.1 )
	&& (R2f>=0.71 && R2f <= 1.5))
	ret= 10;
	if((angA >= 0 && angA <= 60) && (angB >= 60 && angB <=120 )
	&& (angC >= 120 && angC <= 180 ) && (R0f >=0.1 && R0f <= 1.5 )
	&& (R2f>=0.22 && R2f <= 1.5))
	ret= 13;
	if((angA >= 120 && angA <= 180) && (angB >= 0 && angB <=60 )
	&& (angC >= 60 && angC <= 120 ) && (R0f >=0.1 && R0f <= 1.5 )
	&& (R2f>=0.22 && R2f <= 1.5))
	ret= 7;
	if((angA >= 60 && angA <= 120) && (angB >= 120 && angB <=180 )
	&& (angC >= 0 && angC <= 60 ) && (R0f >=0.1 && R0f <= 1.5 )
	&& (R2f>=0.22 && R2f <= 1.5))
	ret= 11;
	if((R0f>=0 && R0f<= 0.1) && (R2f>=0 && R2f <= 0.22))
	ret= 15;
	return ret;
}

int Fault_Selector5(double complex V_pre[], double complex I_pre[], double complex V_pos[], 
					double complex I_pos[], double tol1, double tol2, double R0f, double R2f)
{

	int ret, i, k;
	int var_Z[12];
	double complex Zn_pre[3], Zn_pos[3], Zf_pre[3], Zf_pos[3];
	double complex mZ_pre[6], mZ_pos[6];

	for(i=0;i<12;i++)
	{
		var_Z[i]=2;
	}
	
	ret=0;

	for(i=0; i<3; i++)
	{
		k=i+1;
		Zn_pre[i]=V_pre[i]/I_pre[i];	
		Zn_pos[i]=V_pos[i]/I_pos[i];	

		mZ_pre[i]=Zn_pre[i];
		mZ_pos[i]=Zn_pos[i];
		
		if(i==2)
		{	k=0;}
			
		Zf_pre[i]=(V_pre[i]- V_pre[k]) /(I_pre[i]- I_pre[k]);	
		Zf_pos[i]=(V_pos[i]- V_pos[k]) /(I_pos[i]- I_pos[k]);	
		
		mZ_pre[i+3]=Zf_pre[i];
		mZ_pos[i+3]=Zf_pos[i];
	}

	for(i=0; i<6; i++)
	{
		if (cabs(mZ_pos[i]-mZ_pre[i])>tol1*cabs(mZ_pre[i]))
			{var_Z[i]=1;}
		else
			{var_Z[i]=0;}

		if (cabs(mZ_pos[i]-mZ_pre[i])>tol2*cabs(mZ_pre[i]))
			{var_Z[i+6]=1;}
		else
			{var_Z[i+6]=0;}

	}
	if(R0f+R2f>0.07)
	{
		if( R0f > 0.01)	
		{
		//if(var_Z[9] == 1 && var_Z[8] == 0)
		if(var_Z[9] == 1 && var_Z[8] == 0)
			ret=13;

		//if(var_Z[10] == 1 && var_Z[6] == 0)
		if(var_Z[10] == 1 && var_Z[6] == 0)
			ret=7;

		//if(var_Z[11] == 1 && var_Z[7] == 0)
		if(var_Z[11] == 1 && var_Z[7] == 0)
			ret=11;
	
		if(var_Z[0] == 1 && var_Z[4] == 0)
			ret=9;
	
		if(var_Z[1] == 1 && var_Z[5] == 0)
			ret=5;
	
		if(var_Z[2] == 1 && var_Z[3] == 0)
			ret=3;
		}

		else
		{
		if(var_Z[3] == 1 && var_Z[2] == 0)
			ret=12;

		if(var_Z[4] == 1 && var_Z[0] == 0)
			ret=6;

		if(var_Z[5] == 1 && var_Z[1] == 0)
			ret=10;
		}
	}	
		
		return ret;
}
void PhaseSelectorZ(double complex V_pre[], double complex I_pre[], double complex V_pos[], 
		double complex I_pos[], double tol1, double tol2, double R0f, double R2f, 
		int FS[], double complex dZ[])
{

	int i, k;
	int var_Z[12];
	double complex Zn_pre[3], Zn_pos[3], Zf_pre[3], Zf_pos[3];
	double complex mZ_pre[6], mZ_pos[6];

	for(i=0;i<12;i++)
	{
		var_Z[i]=2;
	}
	for(i=0;i<6;i++)
	{
		FS[i]=0;
	}

	for(i=0; i<3; i++)
	{
		k=i+1;
		Zn_pre[i]=V_pre[i]/I_pre[i];	
		Zn_pos[i]=V_pos[i]/I_pos[i];	

		mZ_pre[i]=Zn_pre[i];
		mZ_pos[i]=Zn_pos[i];
		
		if(i==2)
		{	k=0;}
			
		Zf_pre[i]=(V_pre[i]- V_pre[k]) /(I_pre[i]- I_pre[k]);	
		Zf_pos[i]=(V_pos[i]- V_pos[k]) /(I_pos[i]- I_pos[k]);	
		
		mZ_pre[i+3]=Zf_pre[i];
		mZ_pos[i+3]=Zf_pos[i];
	}
	for(i=0; i<6; i++)
	{
		dZ[i]=mZ_pos[i]-mZ_pre[i];
	}

	for(i=0; i<6; i++)
	{
		if (cabs(mZ_pos[i]-mZ_pre[i])>tol1*cabs(mZ_pre[i]))
			{var_Z[i]=1;}
		else
			{var_Z[i]=0;}

		if (cabs(mZ_pos[i]-mZ_pre[i])>tol2*cabs(mZ_pre[i]))
			{var_Z[i+6]=1;}
		else
			{var_Z[i+6]=0;}

	}
	if(R0f+R2f>0.07)
	{
		if( R0f > 0.01)	
		{
			//if(var_Z[9] == 1 && var_Z[8] == 0)
			if(var_Z[9] == 1 && var_Z[8] == 0)
			{
				FS[3]=1;
			}
			//if(var_Z[10] == 1 && var_Z[6] == 0)
			if(var_Z[10] == 1 && var_Z[6] == 0)
			{	
				FS[4]=1;
			}
			//if(var_Z[11] == 1 && var_Z[7] == 0)
			if(var_Z[11] == 1 && var_Z[7] == 0)
			{
				FS[5]=1;
			}
			if(var_Z[0] == 1 && var_Z[4] == 0)
			{
				FS[0]=1;
				FS[1]=0;
				FS[2]=0;
				FS[3]=0;
				FS[4]=0;
				FS[5]=0;
			}
			if(var_Z[1] == 1 && var_Z[5] == 0)
			{
				FS[0]=0;
				FS[1]=1;
				FS[2]=0;
				FS[3]=0;
				FS[4]=0;
				FS[5]=0;
			}
			if(var_Z[2] == 1 && var_Z[3] == 0)
			{
				FS[0]=0;
				FS[1]=0;
				FS[2]=1;
				FS[3]=0;
				FS[4]=0;
				FS[5]=0;
			}
		}

		else
		{
		if(var_Z[3] == 1 && var_Z[2] == 0)
			FS[3]=1;

		if(var_Z[4] == 1 && var_Z[0] == 0)
			FS[4]=1;

		if(var_Z[5] == 1 && var_Z[1] == 0)
			FS[5]=1;
		}
	}	
		
}

void crearzonasPH(double complex zlinea, double torque, double zo100[], double r0[],
					double x0[], double radio[])
{
	int i;
	double xx1, yy1, a1, b1, c1, d1;
	for(i=0;i<3;i++)
	{
	xx1=cabs(zlinea)/cos(carg(zlinea)-torque*PI/180)*zo100[i]/100*cos(torque/180*PI);
	yy1=cabs(zlinea)/cos(carg(zlinea)-torque*PI/180)*zo100[i]/100*sin(torque/180*PI);
	a1=1;
	b1=-xx1;
	c1=-yy1;
	d1=0;
	r0[i]=-b1/(2*a1);
	x0[i]=-c1/(2*a1);
	radio[i]=sqrt((b1*b1+c1*c1-4*a1*d1)/(4*a1*a1));
	}
}
int zonaPH(double complex zfalla, double r0[], double x0[], double radio[])
{
	int ret=0;
	int i;
	double rfalla=creal(zfalla);
	double xfalla=cimag(zfalla);
	for(i=2;i>=0; i--)
	{
		if((rfalla-r0[i])*(rfalla-r0[i])+(xfalla-x0[i])*(xfalla-x0[i]) <= radio[i]*radio[i])
		ret=i+1;
	}
	return ret;
}
void tiempo(int Nzona, double *timeZ1, double *timeZ2, double *timeZ3, double FrecSamp)
{
	if(Nzona==1)
	{
		*timeZ1=*timeZ1+1/FrecSamp;
		*timeZ2=*timeZ2+1/FrecSamp;
	*timeZ3=*timeZ3+1/FrecSamp;
	}
	if(Nzona==2)
	{
		*timeZ1=0;
		*timeZ2=*timeZ2+1/FrecSamp;
		*timeZ3=*timeZ3+1/FrecSamp;
	}
	if(Nzona==3)
	{
		*timeZ1=0;
		*timeZ2=0;
		*timeZ3=*timeZ3+1/FrecSamp;
	}
	if(Nzona==0)
	{
		*timeZ1=0;
		*timeZ2=0;
		*timeZ3=0;
	}
}
//Carateristica quadrilateral
struct punto
{
	double x,y;
};
struct segmento
{
	struct punto p1;
	struct punto p2;
};
double izq_der(struct segmento S,struct punto P)
{
	double ret;
	struct punto V1;
	struct punto V2;
	V1.x=P.x-S.p1.x;
	V1.y=P.y-S.p1.y;
	V2.x=S.p2.x-S.p1.x;
	V2.y=S.p2.y-S.p1.y;

	ret= V1.x*V2.y-V2.x*V1.y;
	return ret;
}
int isdentro(double X[], double Y[], double x0, double y0, int nP)
{
	int i;
	struct punto P0;
	struct punto P1;
	struct punto P2;
	struct segmento S;
	P0.x=x0;
	P0.y=y0;
	for(i=0;i<nP;i++)
	{
		P1.x=X[i%nP];
		P1.y=Y[i%nP];
		P2.x=X[(i+1)%nP];
		P2.y=Y[(i+1)%nP];
		S.p1=P1;
		S.p2=P2;
		if(izq_der(S,P0)>0)
		return 0;
	}
	return 1;
}
void crearzonasQR(double complex zlinea, double r100[], double x100[],
double XX[][3], double YY[][3], int nVertices[])
{
	int i;
	double R1,X1,RL,XL,cc, L;
	RL=creal(zlinea);
	XL=cimag(zlinea);
	for(i=0;i<3;i++)
	{
		R1=RL*r100[i]/100;
		X1=XL*x100[i]/100;
		L=(X1*tan(9*PI/180)-R1*tan(65*PI/180)*tan(9*PI/180))*sin(81*PI/180)/sin(34*PI/180);
		cc=X1-R1*tan(65*PI/180)-L*cos(25*PI/180);
		XX[0][i]=0;
		XX[1][i]=R1-R1*sin(PI/6)/sin(99*PI/180)*sin(9*PI/180);
		XX[2][i]=R1+X1*tan(9*PI/180);
		nVertices[i]=4;
		if(cc>0)
		{
			XX[3][i]=-R1;
			XX[4][i]=-R1-(X1*tan(9*PI/180)-R1*tan(65*PI/180)*tan(9*PI/180))*sin(81*PI/180)*sin(25*PI/180)/sin(34*PI/180);
			nVertices[i]++;
		}
		else

		{
			XX[3][i]=-X1*tan(25*PI/180);
		}
		YY[0][i]=0;
		YY[1][i]=-R1*sin(PI/6)/sin(99*PI/180)*cos(9*PI/180);
		YY[2][i]=X1;
		YY[3][i]=X1;
		if(cc>0)
		{
			YY[4][i]=(X1*tan(9*PI/180)+R1)*sin(81*PI/180)*cos(25*PI/180)/sin(34*PI/180);
		}
	}
}
int zonaQR(double complex zfalla,double XX[][3],double YY[][3], int nVertices[])
{
	int ret=0;
	int i;
	double rfalla=creal(zfalla);
	double xfalla=cimag(zfalla);
	double X1[5], Y1[5], X2[5], Y2[5], X3[5], Y3[5];
	for(i=0;i<5;i++)
	{
		X1[i]=XX[i][0];
		X2[i]=XX[i][1];
		X3[i]=XX[i][2];
		Y1[i]=YY[i][0];
		Y2[i]=YY[i][1];
		Y3[i]=YY[i][2];
	}
	if(isdentro(X3,Y3,rfalla,xfalla,nVertices[0])==1)
		ret=3;
	if(isdentro(X2,Y2,rfalla,xfalla,nVertices[1])==1)
		ret=2;
	if(isdentro(X1,Y1,rfalla,xfalla,nVertices[2])==1)
		ret=1;
	return ret;
}
#endif
