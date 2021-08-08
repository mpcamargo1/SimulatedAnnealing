#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <unistd.h>
#include <time.h>
#include <math.h>


#define PI 3.14159265359
#define TEMPERATURE 2000
#define MAX_INT 10000
#define MAX_EXT 500
#define BETA 0.978
#define G 15
#define VARIANCE 0.009


/*Utilizado no Ruído Gaussiano*/
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

//Protótipo das funções
double init_temperature();
double init();
double eval(double x);
double nextDouble(long *idum);
double nextGaussian();
double low_temperature(double temperature);
void simulated_annealing(int max_ite_ext,int max_ite_intern,double g);
FILE *fptr;

/*Varíaveis Globais*/

/*Usado no Ruído Gaussiano*/
bool haveNextNextGaussian = false;
double nextNextGuassian;

int main(void ){
	srand(time(0));
	
	fptr = fopen("output.txt","w");
	simulated_annealing(MAX_EXT,MAX_INT,(double) G);

}

double eval(double x){
	
	return pow(2,-2*(pow((x-0.1)/(0.9),2)))*pow(sin(5*PI*x),6);

}
 
double init(){

	return ((double) rand() / (RAND_MAX));

}


double init_temperature(){
	return TEMPERATURE;
}

double low_temperature(double temperature){
	return BETA*temperature;
}


double nextDouble(long *idum){

    int j;
    long k;
    static long iy=0;
    static long iv[NTAB];
    float temp;
 
    if (*idum <= 0 || !iy) {
        if (-(*idum) < 1) *idum=1;
        else *idum = -(*idum);
        for (j=NTAB+7;j>=0;j--) {
            k=(*idum)/IQ;
            *idum=IA*(*idum-k*IQ)-IR*k;
            if (*idum < 0) *idum += IM;
            if (j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }
    k=(*idum)/IQ;
    *idum=IA*(*idum-k*IQ)-IR*k;
    if (*idum < 0) *idum += IM;
    j=iy/NDIV;
    iy=iv[j];
    iv[j] = *idum;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
 }
 

/*Pertubação*/
double nextGaussian(){

	long rand1 = (long) (rand()%RAND_MAX);
	if (haveNextNextGaussian){
		haveNextNextGaussian = false;
		return nextNextGuassian;
	} else {
		double v1,v2,s;
			do {
				v1 = 2*nextDouble(&rand1) - 1;
				rand1 = (long)(rand()%RAND_MAX);
				v2 = 2*nextDouble(&rand1) - 1;
				s = v1*v1 + v2*v2;
			} while (s >= 1 || s == 0);
		nextNextGuassian = v2*sqrt(-2*log(s)/s);
		haveNextNextGaussian = true; 
		return v1*sqrt(-2*log(s)/s);
	}
}


void simulated_annealing(int max_ite_ext,int max_ite_intern,double g){
	double new_x;
	int t_ext,t_int;
	double temperature = init_temperature();
	double x = init();
	printf("Temperatura Inicial: %lf\n",temperature);

	t_ext = 0;

	while(t_ext < max_ite_ext - 1){
		t_int = 0;
		while (t_int < max_ite_intern - 1){
				new_x = x + nextGaussian()*VARIANCE; // Perturba x com ruído guassiano com pequena variância e média zero

				if((new_x >= 0 && new_x <= 1)){ /*Intervalo de interesse [0,1]*/
					if(eval(new_x) > eval(x)) {
						x = new_x;	
						/*Escrever no arquivo*/ /*Somente de 100 em 100 no tempo*/
							if ((MAX_INT*(t_ext)+ t_int)%100==0)fprintf(fptr,"%d %lf\n",MAX_INT*(t_ext)+ t_int,x);
					}
						else {
						if (init(x) < exp( eval(new_x) - eval(x) / (double) temperature)){
							x = new_x;
							/*Escrever no arquivo*/ /*Somente de 100 em 100 no tempo*/
								if ((MAX_INT*(t_ext)+ t_int)%100==0)fprintf(fptr,"%d %lf\n",MAX_INT*(t_ext)+ t_int,x);
						}
					}
				}
				t_int++;
		}
		temperature= low_temperature(temperature);
		t_ext++;
	}
		fclose(fptr);
		printf("Relatório\n Eixo X:%lf\n Temperatura Final:%lf\n",x,temperature);
}
