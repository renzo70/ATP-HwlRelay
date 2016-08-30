#include<math.h>
#include<stdio.h>
#include<complex.h>
#include<time.h>
#include<stdint.h>
#include"./include/relaylib.h"
#include"./include/comtrade.h"

#define PI 3.141592653589793
#define NN 16

//double complex a=cos(2*PI/3)+sin(2*PI/3)*I;

int COMTRADE;
int nA;
int nD;
char nameCOMTRADE[64];
char CTD_station_name[32];
double CTD_IA[NN];
char CTD_time_stamp[NN][32];
char CTD_analog_ID[12][10];
char CTD_digital_ID[12][10];
char local_code[10];
char time_code[10];

struct Analog CTD_vA={.An=1,.ch_id="vA",.ph="a",.ccbm="ccbm",.uu="V",.a=1.0,.b=0,
.skew=10,.min=0,.max=1000,.primary=1000,.secondary=1,.PS='P'};
struct Analog CTD_vB={.An=2,.ch_id="vB",.ph="b",.ccbm="ccbm",.uu="V",.a=1.0,.b=0,
.skew=10,.min=0,.max=1000,.primary=1000,.secondary=1,.PS='P'};
struct Analog CTD_vC={.An=3,.ch_id="vC",.ph="x",.ccbm="ccbm",.uu="V",.a=1.0,.b=0,
.skew=10,.min=0,.max=1000,.primary=1000,.secondary=1,.PS='P'};
struct Digital CTD_FSA={.Dn=1,.ch_id="FSA",.ph="ph",.ccbm="ccbm",.y=0};
struct Digital CTD_FSB={.Dn=1,.ch_id="FSB",.ph="ph",.ccbm="ccbm",.y=0};
struct Digital CTD_FSC={.Dn=1,.ch_id="FSC",.ph="ph",.ccbm="ccbm",.y=0};

struct Analog analog[20];
struct Digital digital[20];

double PT1;
double PT2;
double CT1;
double CT2;
double complex zL0;
double complex zL1;
double lenL;
double SFreq;

int FALTAfull[NN];
int FSfull[NN][6];

double SINPUT[6][NN];

int nsamp[2];

double complex vPHASEpre[3];
double complex iPHASEpre[3];
void c_hwl_relay_m_(double xdata_ar[], double xin_ar[], double xout_ar[], double xvar_ar[])
{
	char TimeNow[30];
	char buffer_time[26];
	char nameCFG[64];
	int it, jt, FALTA; //FALTAfull[4]; 
	int FS[6], FSnow[6];
	double TF, TFF;
	double R0f, R2f;
	double f_ma[6], f_anD[6],f_anR[6] ;
	double complex PHASORS[6]; //Phasor input
	double complex iPHASE[3],iSIME[3], vPHASE[3], vSIME[3];
	double complex ZL0, ZL1;
	double complex dZ[6];
	float ID;
	char ID_str[32];
	FILE *fileCFG;
	double t, starttime, stoptime, fullstep, startstep;

	t=xin_ar[6];


	diffGTM(local_code);
 	diffGTM(time_code);


	ID=xdata_ar[0];
	sprintf(ID_str,"HWL_RELAY-%.0f",ID);
	stoptime=xdata_ar[1];
	fullstep=xdata_ar[2];
	startstep=xdata_ar[3];

	for(it=0;it<6;it++){
		for(jt=NN-1;jt>=1;jt--){
			SINPUT[it][jt]=SINPUT[it][jt-1];
		}
		SINPUT[it][0]=xin_ar[it];
		FSnow[it]=0;
	}
	
	for(it=0;it<NN-1;it++){
		strcpy(CTD_time_stamp[it+1],CTD_time_stamp[it]);
	}
	TimeComtrade(CTD_time_stamp[0]);

//	printf("Muestras leidas\n");
	for(it=0; it<6; it++){
		PHASORS[it]=dft(SINPUT[it], SFreq, 60, 16);
//		rms[i]=RMS(SINPUT[i],16);
		f_ma[it]=cabs(PHASORS[it]);
	}
//	printf("Dft calculada\n");
	for(it=0; it<6; it++)  // Referenciando os sinais a VA
	{
		f_anR[it]=carg(PHASORS[it])-carg(PHASORS[0]);

		if(f_anR[it]<-PI)	f_anR[it]=f_anR[it]+2*PI;
		if(f_anR[it]>PI)	f_anR[it]=f_anR[it]-2*PI;
		
		f_anD[it]=(f_anR[it])*180/PI;
	}

	for(it=0;it<3;it++){
		vPHASE[it]=f_ma[it]*cos(f_anR[it])+f_ma[it]*sin(f_anR[it])*I;
		iPHASE[it]=f_ma[it+3]*cos(f_anR[it+3])+f_ma[it+3]*sin(f_anR[it+3])*I;
	}

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
	
	R0f=fabs(cabs(iSIME[0])/cabs(iSIME[1]));
	R2f=fabs(cabs(iSIME[2])/cabs(iSIME[1]));

	if(R0f+R2f<0.1 || t<1/60.0){
		FALTA=0;
		for(it=0;it<3;it++){
			vPHASEpre[it]=vPHASE[it];
			iPHASEpre[it]=iPHASE[it];
		}
	}
	else{
		FALTA=1;
		for(it=0;it<3;it++){
			vPHASEpre[it]=vPHASEpre[it];
			iPHASEpre[it]=iPHASEpre[it];
		}
		PhaseSelectorZ(vPHASEpre,iPHASEpre,vPHASE,iPHASE,0.01,0.35,R0f,R2f,FSnow,dZ);
	}

//	TF=Fault_Selector5(vPHASEpre,iPHASEpre,vPHASE,iPHASE,0.01,0.35,R0f,R2f);

	for(it=NN-1;it>=1;it--)
	{
		for(jt=0;jt<6;jt++)
		{
			FSfull[it][jt]=FSfull[it-1][jt];
		}
		FALTAfull[it]=FALTAfull[it-1];
	}
	for(jt=0;jt<6;jt++){
		FSfull[0][jt]=FSnow[jt];
	}
	FALTAfull[0]=FALTA;
	
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
	
	//Passando valores

	for(it=0;it<NN;it++){
		CTD_vA.val[it]=SINPUT[0][it];
		CTD_vB.val[it]=SINPUT[1][it];
		CTD_vC.val[it]=SINPUT[2][it];
		CTD_FSA.val[it]=FSfull[it][0];
		CTD_FSB.val[it]=FSfull[it][1];
		CTD_FSC.val[it]=FSfull[it][2];
	}
	analog[0]=CTD_vA;
	analog[1]=CTD_vB;
	analog[2]=CTD_vC;
	digital[0]=CTD_FSA;
	digital[1]=CTD_FSB;
	digital[2]=CTD_FSC;

	if(FALTAfull[1]+FALTAfull[2]+FALTAfull[3]==0)
	{
		if(FALTAfull[0]==1)
		{
			TimeNameComtrade(buffer_time);
			sprintf(nameCOMTRADE,"%s-%s",ID_str,buffer_time);
		
			WriteCFG(nameCOMTRADE,CTD_station_name,"IDE","2013",analog,3,digital,0,
			60,1,1000,20,CTD_time_stamp[NN-1],CTD_time_stamp[0],"ASCII",1,
			time_code,local_code);
			for(it=0;it<NN;it++)
			{
				WriteDAT_init(nameCOMTRADE,analog,3,digital,3);
			}	
		}
	}
	if(FALTAfull[0]+FALTAfull[1]>=2)
	{
		nsamp[0]=nsamp[1]+1;
		WriteDAT(nameCOMTRADE,nsamp[0],analog,3,digital,3);
		nsamp[1]=nsamp[0];
	}

	xout_ar[24]=FS[0];
	xout_ar[25]=FS[1];
//	xout_ar[25]=FALTAfull[0];
	xout_ar[26]=FS[2];
	xout_ar[27]=FS[3];
	xout_ar[28]=FS[4];
	xout_ar[29]=FS[5]; 
	return;
}

void c_hwl_relay_i_(double xdata_ar[], double xin_ar[], double xout_ar[], double xvar_ar[])
{
	int it,jt;
	FILE *Settings;
	double values[10];
	float value;
	int value_int;
	char linea[100];
	char value_str[10];
	char filename[32];
	double R0,X0,R1,X1;
	printf("Lendo o arquivo de configuracao HWL_RELAY-%.0f.set\n",xdata_ar[0]);
	sprintf(filename,"HWL_RELAY-%.0f.set",xdata_ar[0]);
	Settings=fopen(filename,"r");	
	
	printf("Apos ler o arquivo de configuracao HWL_RELAY-%.0f.set\n",xdata_ar[0]);
	fscanf(Settings,"%s", CTD_station_name);
	fgets(linea,100,Settings);

	for(it=0;it<10;it++)
	{	
		fscanf(Settings,"%s", value_str);
		sscanf(value_str,"%f",&value);
		values[it]=value;	
		printf("values[i]: = %f\n",value);
		fgets(linea,100,Settings);
	}
	fscanf(Settings,"%s", value_str);
	sscanf(value_str,"%d",&value_int);
	COMTRADE=value_int;	
	fgets(linea,100,Settings);

	if(COMTRADE==1)
	{
		fscanf(Settings,"%s", value_str);
		sscanf(value_str,"%d",&value_int);
		nA=value_int;	
		fgets(linea,100,Settings);
		
		for(it=0;it<nA;it++)
		{
			fscanf(Settings,"%s", value_str);
			strcpy(CTD_analog_ID[it],value_str);
			printf("%s\n",CTD_analog_ID[it]);
			fgets(linea,100,Settings);
		}

		fscanf(Settings,"%s", value_str);
		sscanf(value_str,"%d",&value_int);
		nD=value_int;	
		fgets(linea,100,Settings);
		
		for(it=0;it<nD;it++)
		{
			fscanf(Settings,"%s", value_str);
			strcpy(CTD_digital_ID[it],value_str);
			fgets(linea,100,Settings);
		}
		
	}	
	
	fclose(Settings);
	
	PT1=values[0];
	PT2=values[1];
	CT1=values[2];
	CT2=values[3];
	X0=values[4];
	R0=values[5];
	R1=values[6];
	X1=values[7];
	lenL=values[8];
	SFreq=values[9];

	zL0=X0+R0*I;	
	zL1=X1+R1*I;	

	nsamp[0]=NN+1;
	nsamp[1]=NN;

	for(jt=0;jt<NN;jt++){
		for(it=0;it<6;it++){
			SINPUT[it][jt]=0;
			FSfull[jt][it]=0;
		}
		FALTAfull[jt]=0;
	}
	for(it=0; it<3; it++){
		vPHASEpre[it]=0+0*I;	
		iPHASEpre[it]=0+0*I;	
	}

	for(it=0; it<30; it++)
		xout_ar[it]=0;

	for(it=0;it<137;it++){
		xvar_ar[it]=0;
	}
//	xvar_ar[114]=1;
//	xvar_ar[115]=1;
//	xvar_ar[116]=1;
	return;
}
