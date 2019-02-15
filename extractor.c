#define MAX(a,b) a>b ? a:b

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <string.h>
#include "trevisan.h"
#include "trevisan.c"


int main(){
	
   	double r;
	double eps;
	double k;
	double alpha;
	int t_req;
	int t;
	int d;
	int n;
	int m;
	int i=0;
	int j=0;
	int l;
	int db;
	double r1;
	size_t **S;
	int x;
	
	char word1[256];
	char word2[256];
	
	FILE *file1;

	FILE *file2;
	
	FILE *file3;
	file3 = fopen("polinomi_irriducibili.txt", "r");

	do{
		printf("\n");
		printf("Enter the file name of the SOURCE (with .txt extension): ");
		scanf("%s", word1);
    	printf("\n");
    	file1 = fopen(word1, "r");
    	
    	if (file1 == NULL){
    	printf("Error: Cannot open input file %s because it couldn't be found or it is empty.\n\n", word1);
    	}
    }while(file1 == NULL);
    
    do{
    	
		printf("Enter the file name of the SEED (with .txt extension): ");
		scanf("%s", word2);
    	printf("\n");
		file2 = fopen(word2, "r");
			
		if (file2 == NULL){
    	printf("Error: Cannot open input file %s because it couldn't be found or it is empty.\n\n\n", word2);
    	}
    }while(file2 == NULL);
		    
	do{
		    
		printf("Enter the min-entropy per bit of source (double): ");
		scanf("%lf", &alpha);	
    	printf("\n");
    	
		if(alpha>1 || alpha<=0){
    		printf("Error: the min-entropy per bit of source must be a number between 0 and 1 (1 can be included).\n\n\n");
    	}
	}while(alpha>1 || alpha<=0);
    
    do{
		printf("Enter the desired error (double): ");
		scanf("%lf", &eps);	
    	printf("\n");
    	
		if(eps>1 || eps<=0){
    		printf("Error: error must be a number between 0 and 1 (1 can be included).\n\n\n");
    	}
	}while(eps>1 || eps<=0);
    
	do{
		printf("Enter the desired weak design you want to use (type in 0 for standard weak design or 1 for block weak design): ");
		scanf("%d", &x);
			
			if(x!=0 && x!=1){
				printf("\n");
				printf("ERROR: enter 0 or 1\n\n\n");
			}
			
	}while(x!=0 && x!=1);
  	
	  printf("\n");
    
	if(x==0){
		r=2*M_E;
	}else{
		r=1;	
	}
	
    fseek(file1, 0L, SEEK_END);
    n = ftell(file1);
    fseek(file1, 0L, SEEK_SET);
    
    
	k=alpha*n;
	m=(int)floor((k-4*log2(1/eps)-6)/r);
	t_req=2*ceil(log2(n)+2*log2(2/eps));
	t=pow(2,ceil(log2(t_req)));                                 
	d=pow(t,2);
	r1=2*M_E;
	l=MAX(1,ceil((log2(m-r1)-log2(t_req-r1))/(log2(r1)-log2(r1-1))));
	db=d*(l+1);
    
    clock_t begin = clock();

	if(m<1){
		printf("Error: with your parameters the random output string bit length is < 1. Change your parameters.\n\n");
		printf("Press enter to continue...");
		fseek(stdin,0,SEEK_END);
		getchar();
		exit(0);
	}
	
	size_t pos = ftell(file1);  
    fseek(file1, 0, SEEK_END);    
    size_t length = ftell(file1); 
    fseek(file1, pos, SEEK_SET);
    
    printf("\n");
    printf("The source length is: %d\n\n", n );
    
	if(length<n){
		printf("Error: the source's bit lenght is longer than the bit length of source input file %s.\n\n", word1);
		printf("Press enter to continue...");
		fseek(stdin,0,SEEK_END);
		getchar();
		exit(0);
	}
	
	printf("\n");
	 
    if(x==0){
    	printf("The seed length required is: %d\n\n", d );
	}else{
		printf("The seed length required is: %d\n\n", db);
	}
	
	printf("\n");
	
	pos = ftell(file2);  
    fseek(file2, 0, SEEK_END);    
    length = ftell(file2);  
	fseek(file2, pos, SEEK_SET);
	
	if (x==0){
		if(length<d){
			printf("Error: the seed's bit lenght is longer than the bit length of seed input file %s.\n\n", word2);
			printf("Press enter to continue...");
			fseek(stdin,0,SEEK_END);
			getchar();
			exit(0);
    	}
	}else{
		if(length<db){
			printf("Error: the seed's bit lenght is longer than the bit length of seed input file %s.\n\n", word2);
			printf("Press enter to continue...");
			fseek(stdin,0,SEEK_END);
			getchar();
			exit(0);
		}
		
	}

	bool *source;
	source=malloc(n*sizeof(bool));

	for (i = 0; i < n; i++){
        fscanf(file1, "%1d", &source[i]);
    }

	bool *seed;
	
	if (x==0){
		
		seed=malloc(d*sizeof(bool));
		
		for (i =0; i <d; i++){
     	   fscanf(file2, "%1d", &seed[i]);
    	}
	}else{
		seed=malloc(db*sizeof(bool));
		
		for (i =0; i <db; i++){
        	fscanf(file2, "%1d", &seed[i]);
    	}
	}
		
	bool **poly_irr;
	
	poly_irr=(bool**)malloc(100*sizeof(bool*));
	for(i=0;i<100;i++){
		poly_irr[i]=malloc(100*sizeof(bool));	
	}
	
	for(i=0; i<100;i++){
		for(j=0;j<100;j++){
			 fscanf(file3, "%1d", &poly_irr[i][j]);
		}
	}
	
	
	S=malloc(m*sizeof(size_t*)); 								
	for(i=0;i<m;i++){
		S[i]=malloc(t_req*sizeof(size_t));	
	}
	
	if(x==0){
		wd(m, t,t_req, S);
	}else{
		bwd(m, t,t_req, S);
	}

	bool *b;
	b=malloc(t_req*sizeof(bool));
	
	bool *rho;
	rho=malloc(m*sizeof(bool));
	
	for(i=0; i<m; i++){
        printf("%.2f%% completed\r",((double)i/(double)m*100));
        fflush(stdout);
		for(j=0; j<t_req; j++){
            
            b[j]=seed[S[i][j]];
		}
		
		rho[i]=one_bit_ext(b,source,n,eps,poly_irr);
		
	}
	
    printf("\n\n");
    
    clock_t end = clock();
    
    char qw[100];
    
    if(x==0){
         strcpy(qw, "random_output_string_WD_");
    }else{
         strcpy(qw, "random_output_string_BWD_");
    }
    
    
    strcat(qw, word1);
    
    FILE *fd;
    fd=fopen(qw, "w");
    
    for(i=0;i<m;i++){
		fprintf(fd, "%d", rho[i]);
	}
    
       
   	printf("The random output string length is: %d\n\n", m);
    
    printf("The random output string is:\n\n");
    
	for(i=0;i<m;i++){
		printf("%d",rho[i]);
	}
    printf("\n\n");
    
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    
    printf("Execution time for the extraction is %.3f seconds.\n\n", time_spent);
    
    char we[100];
    
    if(x==0){
         strcpy(we, "other_data_WD_");
    }else{
         strcpy(we, "other_data_BWD_");
    }
    
    strcat(we, word1);
    
  	FILE *fe;
    fe=fopen(we, "w");
    
    fprintf(fe, "The source length is: %d\n\n", n );
    
    fprintf(fe, "The random output string length is: %d\n\n", m);
    
	if(x==0){
		fprintf(fe, "The seed length required is: %d\n\n", d);
	}else{
		fprintf(fe, "The seed length required is: %d\n\n", db);
	}
	
	fprintf(fe, "Execution time for the extraction is %.3f seconds.\n\n", time_spent);
    
    fclose(fe);
	fclose(fd);
	fclose(file1);
	fclose(file2);
	fclose(file3);
    
	printf("Press enter to continue...");
	fseek(stdin,0,SEEK_END);
	getchar();
	return 0;
}	
