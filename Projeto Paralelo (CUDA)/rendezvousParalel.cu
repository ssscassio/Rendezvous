/**
    CUDA, rendezvousParalel.cu
    Purpose: Calcular as variáveis físicas de um veículo espacial capazes de tornar o Rendezvous possível

    @author Cássio Santos
    @version 1.0 31/07/17 
*/

/**
    Para compilaçao:
        ex.: nvcc rendezvousParalel.cu -o rendezvous
    Para execução:
        ex.: ./rendezvous v-0.005-0.006.dat v-0.007-0.008.dat
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include "cuPrintf.cu"

//Declaração de constantes utilizadas para o calculo dos coeficientes
#define MI 398600.4418
#define EARTH_RADIUS 6378.0
#define N_PRECISION 10

using namespace std;

/**
* Calcular coeficiente A do Rendezvous
* @param N Número de iterações no somatório interno
* @param x0 valor no eixo X da posição relativa inicial entre o satélite e o detrito
* @param y0 valor no eixo Y da posição relativa inicial entre o satélite e o detrito
* @param z0 valor no eixo z da posição relativa inicial entre o satélite e o detrito
* @param xl0 valor no eixo x da velocidade relativa inicial entre o satélite e o detrito
* @param yl0 valor no eixo y da velocidade relativa inicial entre o satélite e o detrito
* @param zl0 valor no eixo z da velocidade relativa inicial entre o satélite e o detrito
* @param Y Gama - Variável física Gama a ser calculado o valor de A
* @param X Chi - Variável física Chi a ser calculado o valor de A
* @param w 
* @param a - 0 - Potência utilizada dentro do somátorio para casos em que o indice do somátorio é utilizado elevado a potências diferentes
*   a = 1  -> n^1
*   a = 2  -> n^2
* @param vex Variável física da Velocidade de exaustão no eixo X a ser calculado o valor de A
* @param vey Variável física da Velocidade de exaustão no eixo Y a ser calculado o valor de A
* @param vez Variável física da Velocidade de exaustão no eixo Z a ser calculado o valor de A
* @returns O coeficiênte A dado os valores iniciais e as variáveis físicas a serem testadas
*/
double __device__ brute_A (int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, double X, double w, int a, double vex, double vey, double vez) {
    double result = 0;
    int n;
    double aux;
    double sum = 0;

    result += (2*xl0)/w - 3*y0 +((2*vex)/w)*log((X+1)/X);
     
    //Calculo do somatorio
    for (n = 1; n <= N; n++) {
        aux = (1/(n*powf(X, n)))*(1/(1+powf(((n*Y)/w),2)))*(((2*vex)/w)+((n*Y*vey)/(w*w)));
        if (n%2 == 0) {
            sum -= aux;
        } else {
            sum += aux;
        }
    }

    result-= sum;
    return result;
}


/**
* Calcular coeficiente B do Rendezvous
* @param N Número de iterações no somatório interno
* @param x0 valor no eixo X da posição relativa inicial entre o satélite e o detrito
* @param y0 valor no eixo Y da posição relativa inicial entre o satélite e o detrito
* @param z0 valor no eixo z da posição relativa inicial entre o satélite e o detrito
* @param xl0 valor no eixo x da velocidade relativa inicial entre o satélite e o detrito
* @param yl0 valor no eixo y da velocidade relativa inicial entre o satélite e o detrito
* @param zl0 valor no eixo z da velocidade relativa inicial entre o satélite e o detrito
* @param Y Gama - Variável física Gama a ser calculado o valor de B
* @param X Chi - Variável física Chi a ser calculado o valor de B
* @param w 
* @param a - 0 - Potência utilizada dentro do somátorio para casos em que o indice do somátorio é utilizado elevado a potências diferentes
*   a = 1  -> n^1
*   a = 2  -> n^2
* @param vex Variável física da Velocidade de exaustão no eixo X a ser calculado o valor de B
* @param vey Variável física da Velocidade de exaustão no eixo Y a ser calculado o valor de B
* @param vez Variável física da Velocidade de exaustão no eixo Z a ser calculado o valor de B
* @returns O coeficiênte B dado os valores iniciais e as variáveis físicas a serem testadas
*/
double __device__ brute_B (int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, double X, double w, int a, double vex, double vey, double vez){
    double result = 0;
    double sum = 0;
    int n;
    double aux;

    result += yl0/w + (vey/w)*log((X+1)/X);

    //Calculo do somatorio
    for (n = 1; n <= N; n++) {
        aux = (1/(n*powf(X,n)))*(1/(1+powf(((n*Y)/w),2)))*(vey/w + (2*n*Y*vex)/(w*w));
        if (n%2 == 0) {//iteração Par
            aux = -aux;
        }
        sum += aux;
    }

    result += sum;
    return result;
}

/**
* Calcular o somatório dos coeficientes Cn do Rendezvous
* @param N Número de iterações no somatório interno
* @param x0 valor no eixo X da posição relativa inicial entre o satélite e o detrito
* @param y0 valor no eixo Y da posição relativa inicial entre o satélite e o detrito
* @param z0 valor no eixo z da posição relativa inicial entre o satélite e o detrito
* @param xl0 valor no eixo x da velocidade relativa inicial entre o satélite e o detrito
* @param yl0 valor no eixo y da velocidade relativa inicial entre o satélite e o detrito
* @param zl0 valor no eixo z da velocidade relativa inicial entre o satélite e o detrito
* @param Y Gama - Variável física Gama a ser calculado o valor de C
* @param X Chi - Variável física Chi a ser calculado o valor de C
* @param w 
* @param a - 0 - Potência utilizada dentro do somátorio para casos em que o indice do somátorio é utilizado elevado a potências diferentes
*   a = 1  -> n^1
*   a = 2  -> n^2
* @param vex Variável física da Velocidade de exaustão no eixo X a ser calculado o valor de C
* @param vey Variável física da Velocidade de exaustão no eixo Y a ser calculado o valor de C
* @param vez Variável física da Velocidade de exaustão no eixo Z a ser calculado o valor de C
* @returns O somatório dos coeficiêntes Cn dado os valores iniciais e as variáveis físicas a serem testadas
*/
double __device__ brute_C (int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, double X, double w, int a, double vex, double vey, double vez){
    double result = 0;
    int n;
    double aux;

    //Calculo do somatorio Cn
    for (n = 1; n <= N; n++) {
        aux = powf(n,a)*vex/(n*powf(X,n)*(1+powf((n*Y/w),2))) + n*Y*powf(n,a)*vey/(n*powf(X,n)*powf(w,2)*(1+powf((n*Y/w),2)));
        if (n%2 == 0) {
            aux = -aux;
        }

        result +=aux;   
    }
    return result;
}

/**
* Calcular coeficiente D do Rendezvous
* @param N Número de iterações no somatório interno
* @param x0 valor no eixo X da posição relativa inicial entre o satélite e o detrito
* @param y0 valor no eixo Y da posição relativa inicial entre o satélite e o detrito
* @param z0 valor no eixo z da posição relativa inicial entre o satélite e o detrito
* @param xl0 valor no eixo x da velocidade relativa inicial entre o satélite e o detrito
* @param yl0 valor no eixo y da velocidade relativa inicial entre o satélite e o detrito
* @param zl0 valor no eixo z da velocidade relativa inicial entre o satélite e o detrito
* @param Y Gama - Variável física Gama a ser calculado o valor de D
* @param X Chi - Variável física Chi a ser calculado o valor de D
* @param w 
* @param a - 0 - Potência utilizada dentro do somátorio para casos em que o indice do somátorio é utilizado elevado a potências diferentes
*   a = 1  -> n^1
*   a = 2  -> n^2
* @param vex Variável física da Velocidade de exaustão no eixo X a ser calculado o valor de D
* @param vey Variável física da Velocidade de exaustão no eixo Y a ser calculado o valor de D
* @param vez Variável física da Velocidade de exaustão no eixo Z a ser calculado o valor de D
* @returns O coeficiênte D dado os valores iniciais e as variáveis físicas a serem testadas
*/
double __device__ brute_D (int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, double X, double w, int a, double vex, double vey, double vez) {
    double result = 0;
    result -= (2*vex* log((X+1)/X))/w;
    result += 4*y0 - 2*xl0/w;
    return result;
}

/**
* Calcular coeficiente E do Rendezvous
* @param N Número de iterações no somatório interno
* @param x0 valor no eixo X da posição relativa inicial entre o satélite e o detrito
* @param y0 valor no eixo Y da posição relativa inicial entre o satélite e o detrito
* @param z0 valor no eixo z da posição relativa inicial entre o satélite e o detrito
* @param xl0 valor no eixo x da velocidade relativa inicial entre o satélite e o detrito
* @param yl0 valor no eixo y da velocidade relativa inicial entre o satélite e o detrito
* @param zl0 valor no eixo z da velocidade relativa inicial entre o satélite e o detrito
* @param Y Gama - Variável física Gama a ser calculado o valor de E
* @param X Chi - Variável física Chi a ser calculado o valor de E
* @param w 
* @param a - 0 - Potência utilizada dentro do somátorio para casos em que o indice do somátorio é utilizado elevado a potências diferentes
*   a = 1  -> n^1
*   a = 2  -> n^2
* @param vex Variável física da Velocidade de exaustão no eixo X a ser calculado o valor de E
* @param vey Variável física da Velocidade de exaustão no eixo Y a ser calculado o valor de E
* @param vez Variável física da Velocidade de exaustão no eixo Z a ser calculado o valor de E
* @returns O coeficiênte E dado os valores iniciais e as variáveis físicas a serem testadas
*/
double __device__ brute_E (int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, double X, double w, int a, double vex, double vey, double vez) {
    double result = 0;

    result -=  3*vex*log((X+1)/X);
    result +=  6*w*y0 - 3*xl0;
    return result;
}

/**
* Calcular o somatório dos coeficientes Fn do Rendezvous
* @param N Número de iterações no somatório interno
* @param x0 valor no eixo X da posição relativa inicial entre o satélite e o detrito
* @param y0 valor no eixo Y da posição relativa inicial entre o satélite e o detrito
* @param z0 valor no eixo z da posição relativa inicial entre o satélite e o detrito
* @param xl0 valor no eixo x da velocidade relativa inicial entre o satélite e o detrito
* @param yl0 valor no eixo y da velocidade relativa inicial entre o satélite e o detrito
* @param zl0 valor no eixo z da velocidade relativa inicial entre o satélite e o detrito
* @param Y Gama - Variável física Gama a ser calculado o valor de Fn
* @param X Chi - Variável física Chi a ser calculado o valor de Fn
* @param w 
* @param a - 0 - Potência utilizada dentro do somátorio para casos em que o indice do somátorio é utilizado elevado a potências diferentes
*   a = 1  -> n^1
*   a = 2  -> n^2
* @param vex Variável física da Velocidade de exaustão no eixo X a ser calculado o valor de Fn
* @param vey Variável física da Velocidade de exaustão no eixo Y a ser calculado o valor de Fn
* @param vez Variável física da Velocidade de exaustão no eixo Z a ser calculado o valor de Fn
* @returns O somatório coeficiênte Fn dado os valores iniciais e as variáveis físicas a serem testadas
*/
double __device__ brute_F(int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, double X, double w, int a, double vex, double vey, double vez) {
    double result = 0;
    double sum = 0;
    int n;
    double aux;


    //Calculo do somatorio
    for (n = 1; n <= N; n++) {
        aux = (1/(n*powf(X,n)))*((2*vey)/w + (4*vex)/(n*Y))/((1+powf((n*Y)/w,2)));
        if (n%2 == 0) {
            aux = - aux;
        }
        aux -= vex/(n*Y);
        aux *= powf(n,a);
        sum += aux;
    }

    result = sum;
    return result;
}

/**
* Calcular coeficiente G do Rendezvous
* @param N Número de iterações no somatório interno
* @param x0 valor no eixo X da posição relativa inicial entre o satélite e o detrito
* @param y0 valor no eixo Y da posição relativa inicial entre o satélite e o detrito
* @param z0 valor no eixo z da posição relativa inicial entre o satélite e o detrito
* @param xl0 valor no eixo x da velocidade relativa inicial entre o satélite e o detrito
* @param yl0 valor no eixo y da velocidade relativa inicial entre o satélite e o detrito
* @param zl0 valor no eixo z da velocidade relativa inicial entre o satélite e o detrito
* @param Y Gama - Variável física Gama a ser calculado o valor de G
* @param X Chi - Variável física Chi a ser calculado o valor de G
* @param w 
* @param a - 0 - Potência utilizada dentro do somátorio para casos em que o indice do somátorio é utilizado elevado a potências diferentes
*   a = 1  -> n^1
*   a = 2  -> n^2
* @param vex Variável física da Velocidade de exaustão no eixo X a ser calculado o valor de G
* @param vey Variável física da Velocidade de exaustão no eixo Y a ser calculado o valor de G
* @param vez Variável física da Velocidade de exaustão no eixo Z a ser calculado o valor de G
* @returns O coeficiênte G dado os valores iniciais e as variáveis físicas a serem testadas
*/
double __device__ brute_G (int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, double X, double w, int a, double vex, double vey, double vez) {
    double result = 0;
    double sum = 0;
    int n;
    double aux;

    result= 2*yl0/w + x0 + (2*vey*(log((X+1)/X)))/w;
    //Calculo do somatorio
    for (n = 1; n <= N; n++) {
        aux = 3*vex/(powf(n,2)*powf(X,n)*w);
        if (n%2 == 0) {
            aux = -aux;
        }
        sum +=aux;
    }
    result-=sum;
    return result;
}

/**
* Calcular coeficiente H do Rendezvous
* @param N Número de iterações no somatório interno
* @param x0 valor no eixo X da posição relativa inicial entre o satélite e o detrito
* @param y0 valor no eixo Y da posição relativa inicial entre o satélite e o detrito
* @param z0 valor no eixo z da posição relativa inicial entre o satélite e o detrito
* @param xl0 valor no eixo x da velocidade relativa inicial entre o satélite e o detrito
* @param yl0 valor no eixo y da velocidade relativa inicial entre o satélite e o detrito
* @param zl0 valor no eixo z da velocidade relativa inicial entre o satélite e o detrito
* @param Y Gama - Variável física Gama a ser calculado o valor de H
* @param X Chi - Variável física Chi a ser calculado o valor de H
* @param w 
* @param a - 0 - Potência utilizada dentro do somátorio para casos em que o indice do somátorio é utilizado elevado a potências diferentes
*   a = 1  -> n^1
*   a = 2  -> n^2
* @param vex Variável física da Velocidade de exaustão no eixo X a ser calculado o valor de H
* @param vey Variável física da Velocidade de exaustão no eixo Y a ser calculado o valor de H
* @param vez Variável física da Velocidade de exaustão no eixo Z a ser calculado o valor de H
* @returns O coeficiênte H dado os valores iniciais e as variáveis físicas a serem testadas
*/
double __device__ brute_H (int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, double X, double w, int a, double vex, double vey, double vez) {
    double result = 0;
    double sum = 0;
    int n;
    double aux;
    
    result = z0;
    //Calculo do somatorio
    for (n = 1; n <= N; n++) {
        aux = ((vez*Y)/(powf(X,n)*powf(w,2)))/(1+powf((n*Y)/w,2));
        if (n%2 == 0) {
            aux = -aux;
        }
        sum += aux;
    }
    result += sum;
    return result;
}

/**
* Calcular coeficiente I do Rendezvous
* @param N Número de iterações no somatório interno
* @param x0 valor no eixo X da posição relativa inicial entre o satélite e o detrito
* @param y0 valor no eixo Y da posição relativa inicial entre o satélite e o detrito
* @param z0 valor no eixo z da posição relativa inicial entre o satélite e o detrito
* @param xl0 valor no eixo x da velocidade relativa inicial entre o satélite e o detrito
* @param yl0 valor no eixo y da velocidade relativa inicial entre o satélite e o detrito
* @param zl0 valor no eixo z da velocidade relativa inicial entre o satélite e o detrito
* @param Y Gama - Variável física Gama a ser calculado o valor de I
* @param X Chi - Variável física Chi a ser calculado o valor de I
* @param w 
* @param a - 0 - Potência utilizada dentro do somátorio para casos em que o indice do somátorio é utilizado elevado a potências diferentes
*   a = 1  -> n^1
*   a = 2  -> n^2
* @param vex Variável física da Velocidade de exaustão no eixo X a ser calculado o valor de I
* @param vey Variável física da Velocidade de exaustão no eixo Y a ser calculado o valor de I
* @param vez Variável física da Velocidade de exaustão no eixo Z a ser calculado o valor de I
* @returns O coeficiênte I dado os valores iniciais e as variáveis físicas a serem testadas
*/
double __device__ brute_I (int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, double X, double w, int a, double vex, double vey, double vez) {
    double result = 0;
    double sum = 0;
    int n;
    double aux;
    
    result = zl0/w - (vez/w)*(log((X+1)/X));

    //Calculo do somatorio
    for (n = 1; n <= N; n++) {
        aux = ((vez)/(powf(n,2)*powf(X,n)*w))/(1+powf((n*Y)/w,2));
        if (n%2 == 0) {
            aux = -aux;
        }
        sum += aux;
    }

    result += sum;
    return result;
}

/**
* Calcular o somatório dos coeficientes Jn do Rendezvous
* @param N Número de iterações no somatório interno
* @param x0 valor no eixo X da posição relativa inicial entre o satélite e o detrito
* @param y0 valor no eixo Y da posição relativa inicial entre o satélite e o detrito
* @param z0 valor no eixo z da posição relativa inicial entre o satélite e o detrito
* @param xl0 valor no eixo x da velocidade relativa inicial entre o satélite e o detrito
* @param yl0 valor no eixo y da velocidade relativa inicial entre o satélite e o detrito
* @param zl0 valor no eixo z da velocidade relativa inicial entre o satélite e o detrito
* @param Y Gama - Variável física Gama a ser calculado o valor de Jn
* @param X Chi - Variável física Chi a ser calculado o valor de Jn
* @param w 
* @param a - 0 - Potência utilizada dentro do somátorio para casos em que o indice do somátorio é utilizado elevado a potências diferentes
*   a = 1  -> n^1
*   a = 2  -> n^2
* @param vex Variável física da Velocidade de exaustão no eixo X a ser calculado o valor de Jn
* @param vey Variável física da Velocidade de exaustão no eixo Y a ser calculado o valor de Jn
* @param vez Variável física da Velocidade de exaustão no eixo Z a ser calculado o valor de Jn
* @returns O somatório coeficiênte Jn dado os valores iniciais e as variáveis físicas a serem testadas
*/
double __device__ brute_J(int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, double X, double w, int a, double vex, double vey, double vez){
    double result = 0;
    double sum = 0;
    int n;
    double aux;

    for (n = 1; n <= N; n++) {
        aux = vez/(n*powf(X,n)*w)/(1+powf((n*Y)/w,2));
        if (n%2 == 0) {
            aux = - aux;
        }
        aux *= powf(n,a);
        sum += aux;
    }
    
    result = sum;
    return result;
}


/**
* Calcular coeficientes A sufixados de 1 a 12 da equação de Rendezvous
* @param A Ponteiro para o array que será modificado e ao fim da execução conterá os valores de A sufixados em cada um de seus indices de 1 à 12
* @param N Número de iterações no somatório interno
* @param x0 valor no eixo X da posição relativa inicial entre o satélite e o detrito
* @param y0 valor no eixo Y da posição relativa inicial entre o satélite e o detrito
* @param z0 valor no eixo z da posição relativa inicial entre o satélite e o detrito
* @param xl0 valor no eixo x da velocidade relativa inicial entre o satélite e o detrito
* @param yl0 valor no eixo y da velocidade relativa inicial entre o satélite e o detrito
* @param zl0 valor no eixo z da velocidade relativa inicial entre o satélite e o detrito
* @param Y Gama - Variável física Gama a ser calculado os coeficientes A sufixados
* @param X Chi - Variável física Chi a ser calculado os coeficientes A sufixados
* @param w 
* @param vex Variável física da Velocidade de exaustão no eixo X a ser calculado os coeficientes A sufixados
* @param vey Variável física da Velocidade de exaustão no eixo Y a ser calculado os coeficientes A sufixados
* @param vez Variável física da Velocidade de exaustão no eixo Z a ser calculado os coeficientes A sufixados
*/
void __device__ brute_all(double *A, int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, int X, double w, double vex, double vey, double vez){

    double a, B, D, e, G, H, I;    
    //Calculando valores de A, B, C, D, E, F, G, H, I e J
    a = brute_A(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 0, vex, vey,vez); 
    B = brute_B(N,x0,y0,z0,xl0,yl0,zl0,  Y, X, w, 0, vex, vey,vez);
    D = brute_D(N,x0,y0,z0,xl0,yl0,zl0,Y, X, w, 0, vex, vey,vez);
    e = brute_E(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 0, vex, vey,vez);
    G = brute_G(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 0, vex, vey,vez);
    H = brute_H(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 0, vex, vey,vez);
    I = brute_I(N,x0,y0,z0,xl0,yl0,zl0,  Y, X, w, 0, vex, vey,vez);
    __syncthreads();

     //Calculando inicialmente a soma (A1² + A3² + A5²)
    A[1]  = 2*a*w + e - Y*brute_F(N,x0,y0,z0, xl0,yl0, zl0, Y, X, w, 1, vex, vey,vez);
    A[3]  = B*w - Y*brute_C(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 1, vex, vey,vez);
    A[5]  = I*w +Y*brute_J(N,x0,y0,z0,xl0, yl0, zl0, Y, X, w, 1, vex, vey,vez);

     //Calculando inicialmente a soma (A2² + A4² + A6²)
    A[2]  = G - 2*B + brute_F(N,x0,y0,z0,xl0,yl0, zl0, Y, X, w, 0, vex, vey,vez);
    A[4]  = a +D + brute_C(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 0, vex, vey,vez);
    A[6]  = H - brute_J(N,x0,y0,z0,xl0,yl0,zl0,  Y, X, w, 0, vex, vey,vez);

     //Calculando inicialmente a soma (A7² + A9² + A11²)
    A[7]  = 2*B*w*w + Y*Y*brute_F(N,x0,y0,z0, xl0, yl0, zl0, Y, X, w, 2, vex, vey,vez);
    A[9]  = Y*Y*brute_C(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 2, vex, vey,vez) - a*w*w;
    A[11] = -w*w*H - Y*Y*brute_J(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 2, vex, vey,vez);

     //Calculando inicialmente a soma (A8² + A10² + A12²) A8=A1 // A10 = A3 // A12 = A5
    A[8]  = A[1];
    A[10] = A[3];
    A[12] = A[5];


}

double __device__ calcularDiferenca (int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double gama, double chi, double w, double vex, double vey, double vez){
    double A[13];
    double a,b,c,d;
    double a1,a2,a3,a4,a5,a6;
    double result;

    brute_all(A, N, x0, y0, z0, xl0, yl0, zl0, gama, chi, w, vex, vey, vez);
    
    //a1 = 1/(A1² + A3² + A5²)
    //a2 = (A2² + A4² + A6²)
    //a3 = 1/(A7² + A9² + A11²)
    //a4 = (A8² + A10² + A12²)
    //a5 = (A1*A2 + A3*A4 + A5*A6)
    //a6 = (A7*A8 + A9*A10 + A11*A12)
    a1 = 1/(A[1]*A[1]+A[3]*A[3]+A[5]*A[5]);
    a2 = A[2]*A[2]+A[4]*A[4]+A[6]*A[6];
    a3 = 1/(A[7]*A[7]+A[9]*A[9]+A[11]*A[11]);
    a4 = A[8]*A[8]+A[10]*A[10]+A[12]*A[12];
    a5 = A[1]*A[2] + A[3]*A[4] + A[5]*A[6];
    a6 = A[7]*A[8] + A[9]*A[10] + A[11]*A[12];

    a = a5*a1;
    b = a1*a2;
    c = a6*a3;
    d = a3*a4;
    
    //Result equivale a diferença entre os dois lados da igualdade que descrevem o Rendezvous
    //Equação do Rendezvous (b-d)² = 4(a-c)(bc-ad)
    result = powf(b-d,2)-4*(a-c)*(b*c-a*d);

    return result;
}

void __global__ calcularRendezvousDevice(double *d_variables){

  double d_x0  = d_variables[0];
  double d_y0  = d_variables[1];
  double d_z0  = d_variables[2];
  double d_xl0 = d_variables[3];
  double d_yl0 = d_variables[4];
  double d_zl0 = d_variables[5];
  double d_w   = d_variables[6];

  double gama = blockIdx.x; //Y(Gama) recebe o valor x atual do bloco;
  double chi  = blockIdx.y; //X(Chi) recebe o valor y atual do bloco;
  double ve   = threadIdx.x; //ve||vex(Velocidade de exaustão) recebe o valor x atual da thread;
  //Conversão de indexs para valores reais
  chi++;
  gama = gama-14;
  gama = powf(10,gama);
  ve++;
  ve = ve/10;

  double yInicial = calcularDiferenca(N_PRECISION, d_x0, d_y0, d_z0, d_xl0, d_yl0, d_zl0, gama, chi, d_w, ve, ve, ve);
  cuPrintf("%.14lf , %lf , %lf , %lf\n",gama, chi, ve, yInicial);
  
}

void calcularRendezvous(double x0, double y0, double z0, double xl0, double yl0, double zl0, double w) {

    int tx = 17; //Número de iterações em Gama (de -14 a 2 com passo 1)
    int ty = 100; //Número de iterações em Chi (de 1 a 100 com passo 1)
    int tz = 100; //Número de iterações em ve (de 0.1 a 10 com passo 0.1)
    dim3 numBlocks(tx,ty);
    size_t size = 7*sizeof(double);
    int threadsPerBlock = tz;


    double *h_variables = (double *)malloc(7*sizeof(double));
    h_variables[0] = x0;
    h_variables[1] = y0;
    h_variables[2] = z0;
    h_variables[3] = xl0;
    h_variables[4] = yl0;
    h_variables[5] = zl0;
    h_variables[6] = w;

    //Alocando a memória das variáveis no Device
    double *d_variables;
    cudaMalloc((void **) &d_variables, size);

    //Copiando as variáveis do host para o Device
    cudaMemcpy(d_variables, h_variables, size, cudaMemcpyHostToDevice);

    calcularRendezvousDevice<<<numBlocks, threadsPerBlock>>>(d_variables);

}

int main(int argc, char **argv){
    //Declaração das variaveis
    double x0,y0,z0; //Componentes das posições relativas
    double xl0,yl0,zl0; //Componentes das velocidades relativas
    double w; //Raio (Constante)
    double raio; //raio (Constante)
    double r0; //Altitude <-----Verificar ----->
    
    //Variaveis inutilizadas (Por enquanto)
    double tempo, alpha, beta, vi, xf, yf, zf, rf, dxf, dyf, vf;

    if(argc == 1){
        printf("Passe o nome dos arquivos de input como parâmetro\n");
        return 1;
    }

    //Aumentando o tamanho do Buffer usado para transferir os dados internos do Device para o Host
    

    // Informação de cuPrintf.cuh // "bufferLen=1048576 1-meg - that's enough for 4096 printfs by all threads put together"
    size_t size = 43520000; //Cada printf necessita de 256 bits 256*100*100*17 
    cudaDeviceSetLimit(cudaLimitPrintfFifoSize, size);

    for(int i = 1; i < argc; i++){ //Tentativa de leitura de cada um dos arquivos passados como parâmetro
        /*Leitura do Arquivo*/
        char * nomeDoArquivo = argv[i];
        FILE *file;
        file = fopen(nomeDoArquivo, "r");
        int b = 0;
        if (file == NULL) { //Verifica se o caminho existe
            break;
        } else {
            while((fscanf(file,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &tempo, &alpha, &beta, &x0, &y0, &z0, &r0, &xl0, &yl0, &zl0, &vi, &xf, &yf, &zf, &rf, &dxf, &dyf, &zf, &vf)) != EOF){ //Enquanto não for o fim do arquivo
                cudaPrintfInit(size);
                
                char nomeDoArquivoDeEscrita[256];
                sprintf( nomeDoArquivoDeEscrita, "%d-output-%d.csv", i, b);
                FILE *fileToWrite;
                fileToWrite = fopen(nomeDoArquivoDeEscrita, "w");

                // Tempo Alpha Beta X0 y0 z0 r0 xl0 yl0 zl0 |Vi| xf yf zf rf dxf dyf dzf |Vf|
                // 456.000000 104 89 -0.725655 2.910444 0.052357 3.000000 0.005108 -0.006719 -0.000104 0.008441 0.000000 0.000000 0.000000 0.000000 -0.001749 -0.005737 -0.000121 0.005999
                raio = EARTH_RADIUS + r0;
                w = sqrt(MI/(raio*raio*raio));
                calcularRendezvous(x0,y0,z0,xl0,yl0,zl0,w);

                // Forçando Flush do buffer do cuPrintf e armazenando em arquivo
                // O segundo parâmetro setado para true faz com que os indices do bloco e da thread que chamou o printf sejam exibidos
                cudaPrintfDisplay(fileToWrite,false);
                cudaPrintfEnd();
		        b++;
            }
        }        
    }

}

