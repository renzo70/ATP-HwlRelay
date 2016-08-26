#include<math.h>
#include<stdio.h>
#include<complex.h>
#include"relaylib.h"
#define PI 3.141592653589793


void c_hwl_relay_m_(double xdata_ar[], double xin_ar[], double xout_ar[], double xvar_ar[])
{
	int it, jt; 
	int FS[6], FSnow[6], FSfull[4][6];
	double TF, TFF;
	double R0f, R2f;
	double f_ma[6], f_anD[6],f_anR[6] ;
	double SINPUT[6][16];
	double complex PHASORS[6]; //Phasor input
	double complex iPHASE[3],iSIME[3], vPHASE[3], vSIME[3], vPHASEpre[3], iPHASEpre[3];
	double complex a=cos(2*PI/3)+sin(2*PI/3)*I;
	double PT1, PT2, CT1, CT2; // Instrument transformers ratios
	double complex zL0, zL1, ZL0, ZL1;
	double lenL, SFreq;
	double complex dZ[6];
	
	PT1=xdata_ar[0];
	PT2=xdata_ar[1];
	CT1=xdata_ar[2];
	CT2=xdata_ar[3];
	zL0=xdata_ar[4]+xdata_ar[5]*I;
	zL1=xdata_ar[6]+xdata_ar[7]*I;
	lenL=xdata_ar[8];
	SFreq=xdata_ar[9];

	//printf("Zonas creadas\n");

    for (it=0; it<6; it++) //guardamos las 6 entradas en el
	{ //arreglo S
		for(jt=0; jt<15; jt++)
		{
			xvar_ar[it*16+jt]=xvar_ar[it*16+jt+1];   // CAMBIAR VARIABLE
		}
		xvar_ar[(it+1)*16-1]=xin_ar[it];
		for(jt=0; jt<16; jt++)
		{
			if(it<3)
				SINPUT[it][jt]=xvar_ar[it*16+jt]*PT1/PT2;
			if(it>=3)
				SINPUT[it][jt]=xvar_ar[it*16+jt]*CT1/CT2;
		}
	}
//	printf("Muestras leidas\n");
	for(it=0; it<6; it++)
	{
		PHASORS[it]=dft(SINPUT[it], SFreq, 60, 16);
//		rms[i]=RMS(SINPUT[i],16);
		f_ma[it]=cabs(PHASORS[it]);
	}
//	printf("Dft calculada\n");
	for(it=0; it<6; it++)  // Referenciando os sinais a VA
	{
		f_anR[it]=carg(PHASORS[it])-carg(PHASORS[0]);
		f_anD[it]=(((carg(PHASORS[it]))-(carg(PHASORS[0]))))*180/PI;
		if(f_anD[it]<-180)
		f_anD[it]=360+f_anD[it];
		if(f_anD[it]>180)
		f_anD[it]=-360+f_anD[it];
	}

	vPHASE[0]=f_ma[0]*cos(f_anR[0])+f_ma[0]*sin(f_anR[0])*I;
	vPHASE[1]=f_ma[1]*cos(f_anR[1])+f_ma[1]*sin(f_anR[1])*I;
	vPHASE[2]=f_ma[2]*cos(f_anR[2])+f_ma[2]*sin(f_anR[2])*I;
	iPHASE[0]=f_ma[3]*cos(f_anR[3])+f_ma[3]*sin(f_anR[3])*I;
	iPHASE[1]=f_ma[4]*cos(f_anR[4])+f_ma[4]*sin(f_anR[4])*I;
	iPHASE[2]=f_ma[5]*cos(f_anR[5])+f_ma[5]*sin(f_anR[5])*I;

	xout_ar[0]=cabs(vPHASE[0]);
	xout_ar[1]=cabs(vPHASE[1]);
	xout_ar[2]=cabs(vPHASE[2]);
	xout_ar[3]=cabs(iPHASE[0]);
	xout_ar[4]=cabs(iPHASE[1]);
	xout_ar[5]=cabs(iPHASE[2]);

	xout_ar[6]=carg(vPHASE[0]);
	xout_ar[7]=carg(vPHASE[1]);
	xout_ar[8]=carg(vPHASE[2]);
	xout_ar[9]=carg(iPHASE[0]);
	xout_ar[10]=carg(iPHASE[1]);
	xout_ar[11]=carg(iPHASE[2]);
	
	simetrica(iPHASE,iSIME);
	simetrica(vPHASE,vSIME);

	vPHASEpre[0]=xvar_ar[96]+xvar_ar[97]*I;
	vPHASEpre[1]=xvar_ar[98]+xvar_ar[99]*I;
	vPHASEpre[2]=xvar_ar[100]+xvar_ar[101]*I;

	iPHASEpre[0]=xvar_ar[102]+xvar_ar[103]*I;
	iPHASEpre[1]=xvar_ar[104]+xvar_ar[105]*I;
	iPHASEpre[2]=xvar_ar[106]+xvar_ar[107]*I;

	R0f=fabs(cabs(iSIME[0])/cabs(iSIME[1]));
	R2f=fabs(cabs(iSIME[2])/cabs(iSIME[1]));

	if(R0f+R2f<0.1)
	{
		xvar_ar[96]=creal(vPHASE[0]);
		xvar_ar[97]=cimag(vPHASE[0]);
		xvar_ar[98]=creal(vPHASE[1]);
		xvar_ar[99]=cimag(vPHASE[1]);
		xvar_ar[100]=creal(vPHASE[2]);
		xvar_ar[101]=cimag(vPHASE[2]);
		xvar_ar[102]=creal(iPHASE[0]);
		xvar_ar[103]=cimag(iPHASE[0]);
		xvar_ar[104]=creal(iPHASE[1]);
		xvar_ar[105]=cimag(iPHASE[1]);
		xvar_ar[106]=creal(iPHASE[2]);
		xvar_ar[107]=cimag(iPHASE[2]);
	}

	TF=Fault_Selector5(vPHASEpre,iPHASEpre,vPHASE,iPHASE,0.01,0.35,R0f,R2f);
	PhaseSelectorZ(vPHASEpre,iPHASEpre,vPHASE,iPHASE,0.01,0.35,R0f,R2f,FSnow,dZ);

	for(it=0;it<18;it++)
	{
		xvar_ar[108+it]=xvar_ar[108+it+6];		
	}
	for(it=0;it<6;it++)
	{
		xvar_ar[126+it]=FSnow[it];	
	}
	for(it=0;it<4;it++)
	{
		for(jt=0;jt<6;jt++)
		{
			FSfull[it][jt]=xvar_ar[108+it*6+jt];
		}
	}

	for(it=0;it<6;it++)
	{
		if(FSfull[0][it]+FSfull[1][it]+ FSfull[2][it]+ FSfull[3][it]==4)
		{
			FS[it]=1;
		}
		else
		{
			FS[it]=0;
		}
	}

	xout_ar[24]=FS[0];
	xout_ar[25]=FS[1];
	xout_ar[26]=FS[2];
	xout_ar[27]=FS[3];
	xout_ar[28]=FS[4];
	xout_ar[29]=FS[5]; 
	return;
}

void c_hwl_relay_i_(double xdata_ar[], double xin_ar[], double xout_ar[], double xvar_ar[])
{
	int it;
	for(it=0; it<30; it++)
		xout_ar[it]=0;

	for(it=0;it<137;it++)
		xvar_ar[it]=0;

	
//	xvar_ar[114]=1;
//	xvar_ar[115]=1;
//	xvar_ar[116]=1;
return;
}

