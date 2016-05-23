#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double ln(double x){
    double result;
    result = log(x)/log(M_E);
}


//1.1 Conhecendo (x0,y0,z0;x'0,y'0,z'0), Y e X calcular ve
void g_ve(float *g_ve, int N, float x0, float y0, float z0, float xl0, float yl0, float zl0, float Y, float X, float w){

    double cvex = 0; //Variavel que armazena o valor do coeficiente de Vex
    double cvey = 0; //Variavel que armazena o valor do coeficiente de Vey
    double cvez = 0; //Variavel que armazena o valor do coeficiente de Vez
    double cont = 0; //Valor livre de variaveis

    // G = cvex*vex + cvey*vey + cvez*vez + cont


//Calculo do somatorio
    int n;
    for(n = 1; n <=N; n++){
        if(n%2 ==0){
            cvex -= 3/(pow(n,2)*pow(X,n)*w);
        }else{
            cvex += 3/(pow(n,2)*pow(X,n)*w);
        }
    }

    cont += 2*yl0/w + x0;
    //Por mudança de base Ln((X+1)/X) = log((X+1)/X)/log(euler)
    cvey += 2*(ln((X+1)/X))/w;

    g_ve[0]= cvex;
    g_ve[1]= cvey;
    g_ve[2]= cvez;
    g_ve[3]= cont;
}


//1.2 Conhecendo (x0,y0,z0;x'0,y'0,z'0), ve e X calcular Y


//1.3 Conhecendo (x0,y0,z0;x'0,y'0,z'0), Y e ve calcular X



int main(){







}