#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

using namespace std;

double __device__ brute_A (int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, int aux2, double w, int a, double vex, double vey, double vez) {
    double result = 0;
    int n;
    double X, aux;
    double sum = 0;
    X = (double) aux2;

    result += (2*xl0)/w - 3*y0 +((2*vex)/w)*log((X+1)/X);
     
    //Calculo do somatorio
    for (n = 1; n <= N; n++) {
        aux = (1/(n*pow(X, n)))*(1/(1+pow(((n*Y)/w),2)))*(((2*vex)/w)+((n*Y*vey)/(w*w)));
        if (n%2 == 0) {
            sum -= aux;
        } else {
            sum += aux;
        }
    }

    result-= sum;
    return result;
}


// Encontrando coeficiente B
double __device__ brute_B (int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, int aux2, double w, int a, double vex, double vey, double vez) {
    double result = 0;
    double sum = 0;
    int n;
    double X,aux;
    X = (double) aux2;

    result += yl0/w + (vey/w)*log((X+1)/X);

    //Calculo do somatorio
    for (n = 1; n <= N; n++) {
        aux = (1/(n*pow(X,n)))*(1/(1+pow(((n*Y)/w),2)))*(vey/w + (2*n*Y*vex)/(w*w));
        if (n%2 == 0) {//iteração Par
            aux = -aux;
        }
        sum += aux;
    }

    result+= sum;

    return result;
}

// Encontrando coeficiente C
double __device__ brute_C (int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, int aux2, double w, int a, double vex, double vey, double vez){
    double result = 0;
    int n;
    double X,aux;
    X = (double) aux2;

    //Calculo do somatorio
    for (n = 1; n <= N; n++) {
        aux = pow(n,a)*vex/(n*pow(X,n)*(1+pow((n*Y/w),2))) + n*Y*pow(n,a)*vey/(n*pow(X,n)*pow(w,2)*(1+pow((n*Y/w),2)));
        if (n%2 == 0) {
            aux = -aux;
        }

        result +=aux;   
    }
    return result;
}

// Encontrando coeficiente D
double __device__ brute_D (int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, int aux2, double w, int a, double vex, double vey, double vez) {
    double result = 0;
    double X;
    X = (double) aux2;

    result -= (2*vex* log((X+1)/X))/w;
    result += 4*y0 - 2*xl0/w;
    return result;
}

// Encontrando coeficiente E
double __device__ brute_E (int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, int aux2, double w, int a, double vex, double vey, double vez) {
    double result = 0;
    double X;
    X = (double) aux2;

    result -=  3*vex*log((X+1)/X);
    result +=  6*w*y0 - 3*xl0;
    return result;
}

// Encontrando coeficiente F
double brute_F(int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, int aux2, double w, int a, double vex, double vey, double vez) {
    double result = 0;
    double sum = 0;
    int n;
    double X,aux;
    X = (double) aux2;


    //Calculo do somatorio
    for (n = 1; n <= N; n++) {
        aux = (1/(n*pow(X,n)))*((2*vey)/w + (4*vex)/(n*Y))/((1+pow((n*Y)/w,2)));
        if (n%2 == 0) {
            aux = - aux;
        }
        aux -= vex/(n*Y);
        aux *= pow(n,a);
        sum += aux;
    }

    result = sum;
    return result;
}

// Encontrando coeficiente G
double __device__ brute_G (int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, int aux2, double w, int a, double vex, double vey, double vez) {
    double result = 0;
    double sum = 0;
    int n;
    double X,aux;
    X = (double) aux2;

    result= 2*yl0/w + x0 + (2*vey*(log((X+1)/X)))/w;
    //Calculo do somatorio
    for (n = 1; n <= N; n++) {
        aux = 3*vex/(pow(n,2)*pow(X,n)*w);
        if (n%2 == 0) {
            aux = -aux;
        }
        sum +=aux;
    }
    result-=sum;
    return result;
}

// Encontrando coeficiente H
double __device__ brute_H (int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, int aux2, double w, int a, double vex, double vey, double vez) {
    double result = 0;
    double sum = 0;
    int n;
    double X,aux;
    X = (double) aux2;
    
    result = z0;
    //Calculo do somatorio
    for (n = 1; n <= N; n++) {
        aux = ((vez*Y)/(pow(X,n)*pow(w,2)))/(1+pow((n*Y)/w,2));
        if (n%2 == 0) {
            aux = -aux;
        }
        sum += aux;
    }
    result += sum;
    return result;
}

// Encontrando coeficiente I
double __device__ brute_I (int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, int aux2, double w, int a, double vex, double vey, double vez) {
    double result = 0;
    double sum = 0;
    int n;
    double X,aux;
    X = (double) aux2;
    
    result = zl0/w - (vez/w)*(log((X+1)/X));

    //Calculo do somatorio
    for (n = 1; n <= N; n++) {
        aux = ((vez)/(pow(n,2)*pow(X,n)*w))/(1+pow((n*Y)/w,2));
        if (n%2 == 0) {
            aux = -aux;
        }
        sum += aux;
    }
    //Por mudança de base log((X+1)/X) = log((X+1)/X)/log(euler)
    result += sum;
    return result;
}

// Encontrando coeficiente J
double __device__ brute_J(int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, int aux2, double w, int a, double vex, double vey, double vez){
    double result = 0;
    double sum = 0;
    int n;
    double X,aux;
    X = (double) aux2;

    for (n = 1; n <= N; n++) {
        aux = vez/(n*pow(X,n)*w)/(1+pow((n*Y)/w,2));
        if (n%2 == 0) {
            aux = - aux;
        }
        aux *= pow(n,a);
        sum += aux;
    }
    
    result = sum;
    return result;
}