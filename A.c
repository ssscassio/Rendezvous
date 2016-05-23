#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double ln(double x){
    double result;
    result = log(x)/log(M_E);
}


//1.1 Conhecendo (x0,y0,z0;x'0,y'0,z'0), Y e X calcular ve
void a_ve(float *a_ve, int N, float x0, float y0, float z0, float xl0, float yl0, float zl0, float Y, float X, float w){

    double cvex = 0; //Variavel que armazena o valor do coeficiente de Vex
    double cvey = 0; //Variavel que armazena o valor do coeficiente de Vey
    double cvez = 0; //Variavel que armazena o valor do coeficiente de Vez
    double cont = 0; //Valor livre de variaveis

    // A = cvex*vex + cvey*vey + cvez*vez + cont


//Calculo do somatorio
    int n;
    for(n = 1; n <=N; n++){
        if(n%2 ==0){
            cvex -= 2/(n*pow(X,n)*w*(1+pow((n*Y/w),2)));

            cvey -= (n*Y)/(n*pow(X,n)*pow(w,2)*(1+pow((n*Y/w),2)));
        }else{
            cvex += 2/(n*pow(X,n)*(1+pow((n*Y/w),2)));

            cvey += (n*Y)/(n*pow(X,n)*pow(w,2)*(1+pow((n*Y/w),2)));
        }
    }

    cont = 2*xl0/w - 3*y0;
    //Por mudança de base Ln((X+1)/X) = log((X+1)/X)/log(euler)
    cvex +=(2*ln((X+1)/X))/w;

    a_ve[0] = cvex;
    a_ve[0] += cvey;
    a_ve[0] += cvez;
    a_ve[1]= cont;

}


//1.2 Conhecendo (x0,y0,z0;x'0,y'0,z'0), ve e X calcular Y
void a_Y(float *a_Y, int N, float x0, float y0, float z0, float xl0, float yl0, float zl0,  float X, float w, float vex, float vey, float vez){

    double cY = 0; //Valor do coeficiente de Y
    double cont = 0; //Valor livre de variaveis

    // A = cvex*vex + cvey*vey + cvez*vez + cont


}


//1.3 Conhecendo (x0,y0,z0;x'0,y'0,z'0), Y e ve calcular X



int main(){







}