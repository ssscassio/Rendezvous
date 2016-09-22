#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double ln(double x){
    double result;
    result = log(x)/log(M_E);
}


//1.1 Conhecendo (x0,y0,z0;x'0,y'0,z'0), Y e X calcular ve
void d_ve(float *d_ve, int N, float x0, float y0, float z0, float xl0, float yl0, float zl0, float Y, float X, float w){

    double cvex = 0; //Variavel que armazena o valor do coeficiente de Vex
    double cvey = 0; //Variavel que armazena o valor do coeficiente de Vey
    double cvez = 0; //Variavel que armazena o valor do coeficiente de Vez
    double cont = 0; //Valor livre de variaveis

    // D = cvex*vex + cvey*vey + cvez*vez + cont


    cvex =  (-2/w)* ln((X+1)/X);

    cont =  4*y0 - 2*xl0/w;

    d_ve[0] = cvex;
    d_ve[0] += cvey;
    d_ve[0] += cvez;
    d_ve[1]= cont;

}


//1.2 Conhecendo (x0,y0,z0;x'0,y'0,z'0), ve e X calcular Y


//1.3 Conhecendo (x0,y0,z0;x'0,y'0,z'0), Y e ve calcular X



int main(){







}