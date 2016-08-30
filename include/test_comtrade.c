#include<stdio.h>
#include"comtrade.h"

int main(int argc, char **argv)
{
	struct Analog analog[1];
	struct Digital digital[1];
	char local_code[7];
	char time_code[7];
	char Time1[30];
	char TimeTrigger[30];

	TimeComtrade(Time1);	
	diffGTM(local_code);
	diffGTM(time_code);
	
	analog[0].An=1;
	strcpy(analog[0].ch_id,"chanel");
	strcpy(analog[0].ph,"ph");
	strcpy(analog[0].ccbm,"ccbm");
	strcpy(analog[0].uu,"kV");
	analog[0].a=3;
	analog[0].b=4;
	analog[0].skew=5;
	analog[0].min=5;
	analog[0].max=5;
	analog[0].primary=1000;
	analog[0].secondary=1;
	analog[0].PS='P';
	
	digital[0].Dn=1;	
	strcpy(digital[0].ch_id,"Chanel_digital");	
	strcpy(digital[0].ph,"ph");	
	strcpy(digital[0].ccbm,"ccbm");	
	digital[0].y=1;	
	
	TimeComtrade(TimeTrigger);	
	WriteCFG(argv[1], "Station_Name","IDE","2013",analog,1,digital,1,60,
	100,1000,20,Time1,TimeTrigger,"ASCII", 1,time_code,local_code);

	return 0;

}

