#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double ln(double x){
    double result;
    result = log(x)/log(M_E);
}


//1.1 Conhecendo (x0,y0,z0;x'0,y'0,z'0), Y e X calcular ve
void c_ve(float *c_ve, int N, float x0, float y0, float z0, float xl0, float yl0, float zl0, float Y, float X, float w, int a){
/**
Usar variavel de entrada para definir se vai ou não multiplicar o N no somatório
pow(n,a)
*/
    double cvex = 0; //Variavel que armazena o valor do coeficiente de Vex
    double cvey = 0; //Variavel que armazena o valor do coeficiente de Vey
    double cvez = 0; //Variavel que armazena o valor do coeficiente de Vez
    double cont = 0; //Valor livre de variaveis

    // Cn = cvex*vex + cvey*vey + cvez*vez + cont


//Calculo do somatorio
    int n;
    for(n = 1; n <=N; n++){
        if(n%2 ==0){
            cvex += pow(n,a)/(n*pow(X,n)*(1+pow((n*Y/w),2)));

            cvey += n*Y*pow(n,a)/(n*pow(X,n)*pow(w,2)*(1+pow((n*Y/w),2)));
        }else{
            cvex -= pow(n,a)/(n*pow(X,n)*(1+pow((n*Y/w),2)));

            cvey -= n*Y*pow(n,a)/(n*pow(X,n)*pow(w,2)*(1+pow((n*Y/w),2)));
        }
    }

    c_ve[0] = cvex;
    c_ve[0] += cvey;
    c_ve[0] += cvez;
    c_ve[1]= cont;

}


//1.2 Conhecendo (x0,y0,z0;x'0,y'0,z'0), ve e X calcular Y


//1.3 Conhecendo (x0,y0,z0;x'0,y'0,z'0), Y e ve calcular X



int main(){







}