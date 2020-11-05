/*
Papanastou Kalliopi
Ergasia: Kataskeuh pssm apo pollaplh stoixish prwteinwn
Input File: fasta format
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//functions
void calcFreq(FILE*);
void calcRelFreq(int);
void calcOdds();
void calcPssm();
void initialize();
void dealloc();
void toPrint(FILE*,double **);
void Print(FILE* ,int **);
void writeToFile(FILE*);

//global variables
char aa[] = {'A','V','L','I','D','E','C','F','G','H','K','N','P','Q','R','S','T','W','Y','M'};
int **freq;
double **relFr,**odds,**pssm;
int length=0;

int main(){
	FILE *fp;
	int n=0,x;
	char c, inputFile[20];
	
	
	printf("Please give the filename: ");
	scanf("%s",&inputFile);
	printf("Press 1 to type intermediate steps else press 0: ");
	scanf("%d", &x);
	
	fp = fopen(inputFile,"r");
	if(fp){		
		printf("File '%s' is opened. \n",inputFile);
	    while ((c = getc(fp)) != EOF){
	    	if(c == '>'){
	    		while((c =getc(fp) )!='\n')continue;
				length = 0;
	    		n++; //# of proteins
	    	}
	    	else if(c != '\n')
				length++; //sequence length
		}    
		fclose(fp);	
	}else{
		printf("File '%s' not found. \n",inputFile);
		return -1;
	}
	printf("Sequence length: %d \nAligned proteins number: %d\n",length,n);
	initialize(length);
	
	fp = fopen(inputFile,"r");
	if(fp){
		printf("File '%s' is opened. \n",inputFile);
		calcFreq(fp);
		fclose(fp);
	}else{
		printf("File '%s' not found. \n",inputFile);
		return -1;
	}
		
	calcRelFreq(n);   
	calcOdds();
	calcPssm();
	
	if(x){
		//frequency results write to file
		fp = fopen ("FrequencyMatrix.txt", "w+");
		if(fp){
			puts("Results in file 'FrequencyMatrix.txt'");
			fprintf(fp,"Protein length: %d \nTotal proteins: %d\n",length,n);
			fprintf(fp,"\t\t\t\t\t\t\t\t\tFREQUENCY MATRIX\n");
			Print(fp,freq);
		}else{
			printf("File 'FrequencyMatrix.txt' not found. \n");
		}
    	fclose(fp);
    	//relative frequency results write to file
		fp = fopen ("RelativeFrequencyMatrix.txt", "w+");
		if(fp){
			puts("Results in file 'RelativeFrequencyMatrix.txt'");
			fprintf(fp,"Protein length: %d \nTotal proteins: %d\n",length,n);
			fprintf(fp,"\n\t\t\t\t\t\t\t\t\tRELATIVE FREQUENCY MATRIX\n");
			toPrint(fp,relFr);
		}else{
			printf("File 'RelativeFrequencyMatrix.txt' not found. \n");
		}
    	fclose(fp);
    	//odds write to file
    	fp = fopen ("OddsMatrix.txt", "w+");
		if(fp){
			puts("Results in file 'OddsMatrix.txt'");
			fprintf(fp,"Protein length: %d \nTotal proteins: %d\n",length,n);
			fprintf(fp,"\n\t\t\t\t\t\t\t\t\tODDS MATRIX\n");
			toPrint(fp,odds);
		}else{
			printf("File 'OddsMatrix.txt' not found. \n");
		}
    	fclose(fp);
    	
	}
	
	
	fp = fopen ("PSSM.txt", "w+");
	
	if(fp){
		puts("Results in file 'PSSM.txt'");
		fprintf(fp,"Protein length: %d \nTotal proteins: %d\n",length,n);
		fprintf(fp,"\n\t\t\t\t\t\t\t\t\tPSSM\n");
		writeToFile(fp);
	}else{
		printf("File 'PSSM.txt' not found. \n");
	}
    fclose(fp);
	
	dealloc();
	
	return 0;
}

void calcFreq(FILE *fp) {
	char c;
	int i,j=0,flag=0;

	while((c = getc(fp)) != EOF){
		if(c=='\n')continue;
		if(c == '>'){
			i=0;
	    	while((c =getc(fp) )!='\n')continue;
		}else{
			for(j=0;j<20;j++){
		 		if(c == aa[j]){
				 	freq[i][j]++;
					flag = 1;
				}						
				if (flag==1){
					flag = 0;
					break;
				}
			}
			i++;
		}	
	}	
}

void calcRelFreq(int n){
	int i,j;
	for(i=0;i<length;i++){
		for(j=0;j<20;j++){
			relFr[i][j] += (double)(freq[i][j])/n;
		}
	}	    
}

void calcOdds(){
	int i,j;
	for(i=0;i<length;i++){
		for(j=0;j<20;j++){
			odds[i][j] = (double)relFr[i][j]/0.05;
		}
	}
}

void calcPssm(){
	int i,j;
	for(i=0;i<length;i++){
		for(j=0;j<20;j++){
			if(odds[i][j] != 0)
				pssm[i][j] = log2(odds[i][j]);
		}
	}
}

void initialize(){
	int i,j;
	freq = (int **)malloc(length*sizeof(*freq));
	relFr = (double **)malloc(length*sizeof(* relFr));
	odds = (double **)malloc(length*sizeof(* odds));
	pssm = (double **)malloc(length*sizeof(* pssm));
	for(j=0;j<length;j++){
    	freq[j] = (int *)malloc(20*sizeof(int));
		relFr[j] = (double *)malloc(20*sizeof(double));
		odds[j] = (double *)malloc(20*sizeof(double));
		pssm[j] = (double *)malloc(20*sizeof(double));
	}

	for(i=0;i<length;i++){
		for(j=0;j<20;j++){
			freq[i][j] = 0;
			relFr[i][j] = 0;
			odds[i][j] = 0;
			pssm[i][j] = 0;
		}
	}
}

void dealloc(){
	int j;
		for(j=0;j<length;j++){
    		free(freq[j]);
			free(relFr[j]);
			free(odds[j]);
			free(pssm[j]);
		}
		free(freq);
		free(relFr);
		free(odds);
		free(pssm);
		
}

void Print(FILE *fp,int ** array){
	int i,j;
	
fprintf(fp,"\t");
	for(i=0;i<20;i++){
		fprintf(fp,"%c\t",aa[i]);
	}
	fprintf(fp,"\n    ================================================================================================================================================================\n");

	for(i=0;i<length;i++){
		if(i<9)
			fprintf(fp,"%4d|\t",i+1);
		else
			fprintf(fp,"%4d|\t",i+1);
		for(j=0;j<20;j++){
			fprintf(fp,"%d\t",array[i][j]);
		}
		fprintf(fp,"\n");
	}
}
	
void toPrint(FILE *fp,double ** array){
	int i,j;
//	double ** ar = (double**)array;
	
	fprintf(fp,"\t");
	for(i=0;i<20;i++){
		fprintf(fp," %c\t",aa[i]);
	}
	fprintf(fp,"\n    ================================================================================================================================================================\n");

	for(i=0;i<length;i++){
		if(i<9)
			fprintf(fp,"%4d|\t",i+1);
		else
			fprintf(fp,"%4d|\t",i+1);
		for(j=0;j<20;j++){
			fprintf(fp,"%.2lf\t",array[i][j]);
		}
		fprintf(fp,"\n");
	}
}

void writeToFile(FILE* fp){
	int i,j;
	fprintf(fp,"\t");
	for(i=0;i<20;i++){
		fprintf(fp," %c\t",aa[i]);
	}
	fprintf(fp,"\n    ================================================================================================================================================================\n");

	for(i=0;i<length;i++){
		if(i<9)
			fprintf(fp,"%4d|\t",i+1);
		else
			fprintf(fp,"%4d|\t",i+1);
		for(j=0;j<20;j++){
			if(pssm[i][j] == 0 )
				fprintf(fp,"-inf\t");
			else
				fprintf(fp,"%.2lf\t",pssm[i][j]);
		}
		fprintf(fp,"\n");
	}
}	
