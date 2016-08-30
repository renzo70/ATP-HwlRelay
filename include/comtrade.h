#ifndef COMTRADE_H
#define COMTRADE_H
#include<math.h>
#include<time.h>
#include<stdio.h>
#include<sys/time.h>
#include<string.h>
#include<stdbool.h>

#define NN 16

struct Analog{
	int An;
	char ch_id[128];
	char ph[2];
	char  ccbm[64];
	char  uu[32];
	double a;
	double b;
	double  skew;
	double  min;
	double max;
	double  primary;
	double  secondary;
	char PS;
	double val[NN];
};
struct Digital{
	int Dn; 
	char ch_id[128];
	char  ph[2];
	char ccbm[64];
	bool y;
	bool val[NN];
};

void TimeComtrade(char Time[])
{
	time_t rawtime;
	struct tm *timeinfo;
	struct timeval tv;
	char buffer[26];
	int microsec;

	gettimeofday(&tv,NULL);
	
	timeinfo=localtime(&tv.tv_sec);
	strftime(buffer, 26,"%d/%m/%Y,%H:%M:%S",timeinfo);

	sprintf(Time,"%s.%06d",buffer,tv.tv_usec);
}
void TimeNameComtrade(char Time[])
{
	time_t rawtime;
	struct tm *timeinfo;
	struct timeval tv;
	int microsec;

	gettimeofday(&tv,NULL);
	
	timeinfo=localtime(&tv.tv_sec);
	strftime(Time, 26,"%Y%m%d%H%M%S",timeinfo);
}

void diffGTM(char diff[])
{
	time_t rawtime;
   	struct tm iGMT;
   	struct tm iLocal;
	
	char dH_str[4];
	char dM_str[3];
	int dMT,dH,dM,dMabs;

   	rawtime=time(NULL);
   /* Get GMT time */
   	gmtime_r(&rawtime,&iGMT);
   	localtime_r(&rawtime,&iLocal);

	dMT=iLocal.tm_min-iGMT.tm_min+(iLocal.tm_hour - iGMT.tm_hour)*60+
		(iLocal.tm_mday-iGMT.tm_mday)*60*24+
		(iLocal.tm_mon-iGMT.tm_mon)*60*24*30;
	dMabs=abs(dMT);
	dM=dMabs%60;
	dH=dMabs/60;
	
	if(dMT<0)
	{
		if(dM!=0)
		{
			sprintf(diff,"-%dh%02d",dH,dM);
		}
		else
		{
			sprintf(diff,"-%d",dH);
		}
	}
	else
	{
		if(dM!=0)
		{
			sprintf(diff,"%dh%02d",dH,dM);
		}
		else
		{
			sprintf(diff,"%d",dH);
		}
	}

}

void WriteCFG(char filename[],
 char station_name[], char rec_dev_id[], char rev_year[],
struct Analog analog[],int nA, struct Digital digital[], int nD,
double If, int nrates, double samp, int endsamp,
char Time1[], char TimeTrigger[],char ft[], double timemult,
char time_code[], char local_code[])
{
	char Time[40];
	char filenamefull[200];
	FILE *fileCFG;
	int it, jt;
	int TT;
	TT=nA+nD;

	sprintf(filenamefull,"%s.cfg",filename);
	fileCFG=fopen(filenamefull,"w");

	fprintf(fileCFG, "%s,%s,%s\n", station_name, rec_dev_id, rev_year);
	fprintf(fileCFG, "%d,%d,%d\n", TT,nA,nD);
	
	for(it=0;it<nA;it++)
	{
		fprintf(fileCFG, "%d,%s,%s,%s,%s,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%c\n",
		analog[it].An,analog[it].ch_id,analog[it].ph,analog[it].ccbm,analog[it].uu,
		analog[it].a,analog[it].b,analog[it].skew, analog[it].min, analog[it].max, 
		analog[it].primary,analog[it].secondary,analog[it].PS);
	}
	for(it=0;it<nD;it++)
	{
		fprintf(fileCFG, "%d,%s,%s,%s,%d\n",
		digital[it].Dn, digital[it].ch_id, digital[it].ph, digital[it].ccbm,digital[it].y);
	}
	fprintf(fileCFG,"%.2f\n",If); //if; I uppercase
	fprintf(fileCFG,"%d\n",nrates);
	fprintf(fileCFG,"%.2f,%d\n",samp,endsamp);

	fprintf(fileCFG,"%s\n",Time1);
	fprintf(fileCFG,"%s\n",TimeTrigger);
	fprintf(fileCFG,"%s\n",ft);
	fprintf(fileCFG,"%f\n",timemult);
	fprintf(fileCFG,"%s,%s\n",time_code,local_code);
	fprintf(fileCFG,"9,3\n");
	fclose(fileCFG);
}

void WriteDAT_init(char filename[],struct Analog analog[], int nA, 
				struct Digital digital[], int nD)
{
	int it,jt;
	char filenamefull[200];
	FILE *fileDAT;
	sprintf(filenamefull,"%s.dat",filename);
	fileDAT=fopen(filenamefull,"w");

	for(jt=0;jt<NN;jt++){
		fprintf(fileDAT,"%d,",jt+1);
	
		for(it=0;it<nA;it++){
			fprintf(fileDAT,",%f",analog[it].val[NN-1-jt]);
		}
		for(it=0;it<nD;it++){
			fprintf(fileDAT,",%d",digital[it].val[NN-1-jt]);
		}
		fprintf(fileDAT,"\n");
	}

	fclose(fileDAT);
}
void WriteDAT(char filename[],int sample_number,
				struct Analog analog[], int nA, 
				struct Digital digital[],int nD)
{
	int it,jt;
	char filenamefull[200];
	FILE *fileDAT;
	sprintf(filenamefull,"%s.dat",filename);
	fileDAT=fopen(filenamefull,"a");
	
	fprintf(fileDAT,"%d,", sample_number);
	
	for(it=0;it<nA;it++){
		fprintf(fileDAT,",%f",analog[it].val[0]);
	}
	for(it=0;it<nD;it++){
		fprintf(fileDAT,",%d",digital[it].val[0]);
	}
	fprintf(fileDAT,"\n");
	fclose(fileDAT);
}
#endif
