#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double ln(double x){
    double result;
    result = log(x)/log(M_E);
}


//1.1 Conhecendo (x0,y0,z0;x'0,y'0,z'0), Y e X calcular ve
void i_ve(float *i_v, int N, float x0, float y0, float z0, float xl0, float yl0, float zl0, float Y, float X, float w){

    double cvex = 0; //Variavel que armazena o valor do coeficiente de Vex
    double cvey = 0; //Variavel que armazena o valor do coeficiente de Vey
    double cvez = 0; //Variavel que armazena o valor do coeficiente de Vez
    double cont = 0; //Valor livre de variaveis

    // B = cvex*vex + cvey*vey + cvez*vez + cont


//Calculo do somatorio
    int n;
    for(n = 1; n <=N; n++){
        if(n%2 ==0){
            cvez += 1/(pow(n,2)*pow(X,n)*w*(1+pow((n*Y/w),2)));
        }else{
            cvez -= 1/(pow(n,2)*pow(X,n)*w*(1+pow((n*Y/w),2)));
        }
    }

    cont += zl0/w;
    //Por mudança de base Ln((X+1)/X) = log((X+1)/X)/log(euler)
    cvez -= (ln((X+1)/X))/w;

    //vex = vez = vey
    //ve = sqrt(3v*v)
    i_ve[0] = cvex;
    i_ve[0] += cvey;
    i_ve[0] += cvez;
    i_ve[1]= cont;
}


//1.2 Conhecendo (x0,y0,z0;x'0,y'0,z'0), ve e X calcular Y


//1.3 Conhecendo (x0,y0,z0;x'0,y'0,z'0), Y e ve calcular X



int main(){







}