#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double ln(double x){
    double result;
    result = log(x)/log(M_E);
}


//1.1 Conhecendo (x0,y0,z0;x'0,y'0,z'0), Y e X calcular ve
void f_ve(float *f_ve, int N, float x0, float y0, float z0, float xl0, float yl0, float zl0, float Y, float X, float w, int a){
/**
Usar variavel de entrada para definir se vai ou não multiplicar o N no somatório
pow(n,a)
*/
    double cvex = 0; //Variavel que armazena o valor do coeficiente de Vex
    double cvey = 0; //Variavel que armazena o valor do coeficiente de Vey
    double cvez = 0; //Variavel que armazena o valor do coeficiente de Vez
    double cont = 0; //Valor livre de variaveis

    // Fn = cvex*vex + cvey*vey + cvez*vez + cont


//Calculo do somatorio
    int n;
    for(n = 1; n <=N; n++){
        if(n%2 ==0){
            cvex += pow(n,a)*4/(n*pow(X,n)*n*Y*(1+pow((n*Y/w),2)));

            cvey += 2*pow(n,a)/(n*pow(X,n)*w*(1+pow((n*Y/w),2)));

        }else{
            cvex -= pow(n,a)*4/(n*pow(X,n)*n*Y*(1+pow((n*Y/w),2)));

            cvey -= 2*pow(n,a)/(n*pow(X,n)*w*(1+pow((n*Y/w),2)));
        }
    }

    cvex -= 1/n*Y;

    f_ve[0] = cvex;
    f_ve[0] += cvey;
    f_ve[0] += cvez;
    f_ve[1]= cont;

}


//1.2 Conhecendo (x0,y0,z0;x'0,y'0,z'0), ve e X calcular Y


//1.3 Conhecendo (x0,y0,z0;x'0,y'0,z'0), Y e ve calcular X



int main(){







}