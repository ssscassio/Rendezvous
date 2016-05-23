#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double ln(double x){
    double result;
    result = log(x)/log(M_E);
}


//1.1 Conhecendo (x0,y0,z0;x'0,y'0,z'0), Y e X calcular ve
void e_ve(float *e_ve, int N, float x0, float y0, float z0, float xl0, float yl0, float zl0, float Y, float X, float w){

    double cvex = 0; //Variavel que armazena o valor do coeficiente de Vex
    double cvey = 0; //Variavel que armazena o valor do coeficiente de Vey
    double cvez = 0; //Variavel que armazena o valor do coeficiente de Vez
    double cont = 0; //Valor livre de variaveis

    // E = cvex*vex + cvey*vey + cvez*vez + cont


    cvex =  -3* ln((X+1)/X);

    cont =  6*w*y0 - 3*xl0;

    e_ve[0] = cvex;
    e_ve[0] += cvey;
    e_ve[0] += cvez;
    e_ve[1]= cont;

}


//1.2 Conhecendo (x0,y0,z0;x'0,y'0,z'0), ve e X calcular Y


//1.3 Conhecendo (x0,y0,z0;x'0,y'0,z'0), Y e ve calcular X



int main(){







}