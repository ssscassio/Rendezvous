/**
    Para compilaçao:
        ex.: gcc -o rendezvous RendezvousZeroFinding.c -lm -std=c99
    Para execução:
        ex.: ./rendezvous v-0.005-0.006.dat v-0.007-0.008.dat >saida.dat
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//Declaração de constantes
#define MI 398600.4418
#define EARTH_RADIUS 6378.0
#define _USE_MATH_DEFINES
#define E 2.71828
//Funções usadas para calculo de Rendezvous
double ln (double x) {
    double result;

    result = log(x)/log(E);

    return result;
}

//Encontrando coeficiente A
double brute_A (int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, int X, double w, int a, double vex, double vey, double vez) {
    double result = 0;
    int n;

    //Calculo do somatorio
    for (n = 1; n <= N; n++) {
        if (n%2 == 0) {
            result -= 2*vex/(n*pow(X,n)*w*(1+pow((n*Y/w),2)));
            result -= (n*Y*vey)/(n*pow(X,n)*pow(w,2)*(1+pow((n*Y/w),2)));
        } else {
            result += 2*vex/(n*pow(X,n)*(1+pow((n*Y/w),2)));
            result += (n*Y*vey)/(n*pow(X,n)*pow(w,2)*(1+pow((n*Y/w),2)));
        }
    }
    result += 2*xl0/w - 3*y0;
    //Por mudança de base Ln((X+1)/X) = log((X+1)/X)/log(euler)
    result += (2*vex*ln((X+1)/X))/w;

    return result;
}

// Encontrando coeficiente B
double brute_B (int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, int X, double w, int a, double vex, double vey, double vez) {
    double result = 0;
    int n;

    //Calculo do somatorio
    for (n = 1; n <= N; n++) {
        if (n%2 == 0) {//iteração Par
            result += 2*vex*n*Y/(n*pow(X,n)*pow(w,2)*(1+pow((n*Y/w),2)));
            result += vey/(n*pow(X,n)*w*(1+pow((n*Y/w),2)));
        } else {//iteração Impar
            result -= 2*n*vex*Y/(n*pow(X,n)*pow(w,2)*(1+pow((n*Y/w),2)));
            result -= vey/(n*pow(X,n)*w*(1+pow((n*Y/w),2)));
        }
    }
    result += yl0/w;
    //Por mudança de base Ln((X+1)/X) = log((X+1)/X)/log(euler)
    result += vey*(ln((X+1)/X))/w;

    return result;
}

// Encontrando coeficiente C
double brute_C (int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, int X, double w, int a, double vex, double vey, double vez){
    double result = 0;
    int n;

    //Calculo do somatorio
    for (n = 1; n <= N; n++) {
        if (n%2 == 0) {
            result += pow(n,a)*vex/(n*pow(X,n)*(1+pow((n*Y/w),2)));
            result += n*Y*pow(n,a)*vey/(n*pow(X,n)*pow(w,2)*(1+pow((n*Y/w),2)));
        } else {
            result -= pow(n,a)*vex/(n*pow(X,n)*(1+pow((n*Y/w),2)));
            result -= n*Y*pow(n,a)*vey/(n*pow(X,n)*pow(w,2)*(1+pow((n*Y/w),2)));
        }
    }

    return result;
}

// Encontrando coeficiente D
double brute_D (int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, int X, double w, int a, double vex, double vey, double vez) {
    double result = 0;

    result -= (2*vex* ln((X+1)/X))/w;
    result += 4*y0 - 2*xl0/w;

    return result;
}

// Encontrando coeficiente E
double brute_E (int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, int X, double w, int a, double vex, double vey, double vez) {
    double result = 0;

    result -=  3*vex*ln((X+1)/X);
    result +=  6*w*y0 - 3*xl0;

    return result;
}

// Encontrando coeficiente F
double brute_F(int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, int X, double w, int a, double vex, double vey, double vez) {
    double result = 0;
    int n;

    //Calculo do somatorio
    for (n = 1; n <= N; n++) {
        if (n%2 ==0) {
            result += pow(n,a)*vex*4/(n*pow(X,n)*n*Y*(1+pow((n*Y/w),2)));
            result += 2*vey*pow(n,a)/(n*pow(X,n)*w*(1+pow((n*Y/w),2)));
        } else {
            result -= pow(n,a)*vex*4/(n*pow(X,n)*n*Y*(1+pow((n*Y/w),2)));
            result -= 2*vey*pow(n,a)/(n*pow(X,n)*w*(1+pow((n*Y/w),2)));
        }
    }
    result -= vex/n*Y;

    return result;
}

// Encontrando coeficiente G
double brute_G (int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, int X, double w, int a, double vex, double vey, double vez) {
    double result = 0;
    int n;

    //Calculo do somatorio
    for (n = 1; n <= N; n++) {
        if (n%2 == 0) {
            result -= 3*vex/(pow(n,2)*pow(X,n)*w);
        } else {
            result += 3*vex/(pow(n,2)*pow(X,n)*w);
        }
    }
    result+= 2*yl0/w + x0;
    //Por mudança de base Ln((X+1)/X) = log((X+1)/X)/log(euler)
    result+= 2*vey*(ln((X+1)/X))/w;

    return result;
}

// Encontrando coeficiente H
double brute_H (int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, int X, double w, int a, double vex, double vey, double vez) {
    double result = 0;
    int n;
    
    //Calculo do somatorio
    for (n = 1; n <= N; n++) {
        if (n%2 == 0) {
            result += vez*Y/(pow(w,2)*pow(X,n)*(1+pow((n*Y/w),2)));
        } else {
            result -= vez*Y/(pow(w,2)*pow(X,n)*(1+pow((n*Y/w),2)));
        }
    }
    result += z0;

    return result;
}

// Encontrando coeficiente I
double brute_I (int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, int X, double w, int a, double vex, double vey, double vez) {
    double result = 0;
    int n;
    
    //Calculo do somatorio
    for (n = 1; n <= N; n++) {
        if (n%2 == 0) {
            result += vez/(pow(n,2)*pow(X,n)*w*(1+pow((n*Y/w),2)));
        } else {
            result -= vez/(pow(n,2)*pow(X,n)*w*(1+pow((n*Y/w),2)));
        }
    }
    result += zl0/w;
    //Por mudança de base Ln((X+1)/X) = log((X+1)/X)/log(euler)
    result -= vez*(ln((X+1)/X))/w;

    return result;
}

// Encontrando coeficiente J
double brute_J(int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, int X, double w, int a, double vex, double vey, double vez){
    double result = 0;
    int n;

    //Calculo do somatorio
    for(n = 1; n <= N; n++){
        if (n%2 == 0) {
            result += vez*pow(n,a)/(n*pow(X,n)*(1+pow((n*Y/w),2)));
        } else {
            result -= vez*pow(n,a)/(n*pow(X,n)*(1+pow((n*Y/w),2)));
        }
    }

    return result;
}

double brute_all(double *A, int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, int X, double w, double vex, double vey, double vez){


    double result = 0;

     //Calculando inicialmente a soma (A1² + A3² + A5²)
    A[1]  = 2*brute_A(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 0, vex, vey,vez)*w + brute_E(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 0, vex, vey,vez) - Y*brute_F(N,x0,y0,z0, xl0,yl0, zl0, Y, X, w, 1, vex, vey,vez);
    A[3]  = brute_B(N,x0,y0,z0,xl0,yl0,zl0,  Y, X, w, 0, vex, vey,vez)*w  - Y*brute_F(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 0, vex, vey,vez);
    A[5]  = brute_I(N,x0,y0,z0,xl0,yl0,zl0,  Y, X, w, 0, vex, vey,vez)*w  - Y*brute_F(N,x0,y0,z0,xl0,yl0,zl0,  Y, X, w, 0, vex, vey,vez)*w+Y*brute_J(N,x0,y0,z0,xl0, yl0, zl0, Y, X, w, 0, vex, vey,vez)*w  - Y*brute_F(N,x0,y0,z0,xl0,yl0, zl0, Y, X, w, 1, vex, vey,vez);;

     //Calculando inicialmente a soma (A2² + A4² + A6²)
    A[2]  = brute_G(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 0, vex, vey,vez) - 2*brute_B(N,x0,y0,z0,xl0,yl0,zl0,  Y, X, w, 0, vex, vey,vez) + brute_F(N,x0,y0,z0,xl0,yl0, zl0, Y, X, w, 0, vex, vey,vez);
    A[4]  = brute_A(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 0, vex, vey,vez) +brute_D(N,x0,y0,z0,xl0,yl0,zl0,Y, X, w, 0, vex, vey,vez) + brute_C(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 0, vex, vey,vez);
    A[6]  = brute_H(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 0, vex, vey,vez) - brute_J(N,x0,y0,z0,xl0,yl0,zl0,  Y, X, w, 0, vex, vey,vez);

     //Calculando inicialmente a soma (A7² + A9² + A11²)
    A[7]  = 2*brute_B(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 0, vex, vey,vez)*w*w + Y*Y*brute_F(N,x0,y0,z0, xl0, yl0, zl0, Y, X, w, 2, vex, vey,vez);
    A[9]  = Y*Y*brute_C(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 2, vex, vey,vez) - brute_A(N,x0,y0,z0,xl0,yl0,zl0,  Y, X, w, 0, vex, vey,vez)*w*w;
    A[11] = -w*w*brute_H(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 0, vex, vey,vez) - Y*Y*brute_J(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 2, vex, vey,vez);

     //Calculando inicialmente a soma (A8² + A10² + A12²) A8=A1 // A10 = A3 // A12 = A5
    A[8]  = A[1];
    A[10] = A[3];
    A[12] = A[5];
}

double calcularDiferenca (int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double gama, double chi, double w, double vex, double vey, double vez){
    double A[13];
    double a,b,c,d;
    double a1,a2,a3,a4,a5,a6;

    brute_all(A, 10, x0, y0, z0, xl0, yl0, zl0, gama, chi, w, vex, vey, vez);
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

    return pow(b-d,2)-4*(a-c)*(b*c-a*d);
}

void calcularRendezvousBisseccao(double x0, double y0, double z0, double xl0, double yl0, double zl0, double w) {
    double vex,vey,vez; //Variaveis Iteraveis //Componentes das Velocidades de exaustao
    double ve; //Modulo da velocidade de exaustao
    double chi; //Variavel Iteravel
    double gama; //Variavel Iteravel

    //Variaveis de iteração
    int i;
    int iteracoes;

    //Variaveis para iteração por bissecção
    //X diz respeito ao valor de Ve sendo Iterado
    //Y diz respeito ao resultado da diferença
    double xInicial, xFinal, xMedio, yMedio, yInicial, yFinal;

    /*Parte para leitura das variaveis nos arquivo<------A fazer ----->*/
    printf("Valores de entrara: x0 y0 z0 xl0 yl0 zl0 w\n");
    printf("%lf %lf %lf %lf %lf %lf\n\n", x0,y0,z0, xl0, yl0, zl0);

    for (i = -14; i <= 2 ; i++) {//Iterando gama
        gama = pow(10,i);
        for(chi = 1; chi <= 100 ; chi++) {//Iterando Chi
            for(vex = 0.1; vex <= 10 ; vex+=0.1) {//Iterando nas componentes das velocidades de exaustão
                xInicial = vex;
                xFinal = vex+0.1;
                
                //Calculando diferença no ponto xInicial
                vex = vey = vez = xInicial;
                yInicial = calcularDiferenca(10, x0, y0, z0, xl0, yl0, zl0, gama, chi, w, vex, vey, vez);

                //Calculando diferença no ponto xFinal
                vex = vey = vez = xFinal;
                yFinal = calcularDiferenca(10, x0, y0, z0, xl0, yl0, zl0, gama, chi, w, vex, vey, vez);

                //Verifica se existe uma mudança de sinal entre os dois pontos 
                //(Indicativo de possivel raiz)
                if(yInicial*yFinal <=0) {
                    iteracoes = 0;
                    
                    do { //Iterar até o limite de iterações ou até achar um valor exatamete com precisão 10^-4
                        xMedio = (xInicial+xFinal)/2; //Pegar o valor de X entre os dois limites
                        vex = vey = vez = xMedio;
                        yMedio = calcularDiferenca(10, x0, y0, z0, xl0, yl0, zl0, gama, chi, w, vex, vey, vez);
                        if( yMedio*yInicial < 0) {//Mudança de sinal ocorre entre o limite inferior e o ponto medio
                            xFinal = xMedio;
                            yFinal = yMedio;
                        } else { //Mudança de sinal ocorre entre o limite superior e o ponto medio
                            xInicial = xMedio;
                            yInicial = yMedio;
                        }
                        iteracoes++;
                    } while(abs(yMedio) > 0.0001 && iteracoes <20);
                    //Conjuntos de valoes no final do laço dizem respeito a um conjunto ideal
                    printf("Encontrou yMedio: %lf, gama: %lf, vex: %lf, chi %lf \n", yMedio, gama, xMedio, chi); //<----Mudar por salvar valores em arquivo---->
                } else if (abs(yInicial) < 0.0001) {
                    printf("Encontrou yInicial: %lf, gama: %lf, vex: %lf, chi %lf \n", yInicial, gama, xInicial, chi); //<----Mudar por salvar valores em arquivo---->
                }
            }
        }
    }
}

int main (int argc, char **argv){
    //Declaração das variaveis
    double x0,y0,z0; //Componentes das posições relativas
    double xl0,yl0,zl0; //Componentes das velocidades relativas
    double w; //Raio (Constante)
    double raio; //raio (Constante)
    double r0; //Altitude <-----Verificar ----->
    
    //Variaveis inutilizadas (Por enquanto)
    double tempo, alpha, beta, vi, xf, yf, zf, rf, dxf, dyf, dzf, vf;

    if(argc == 1){
        printf("Passe o nome dos arquivos de input como parâmetro\n");
        return 1;
    }

    for(int i = 1; i < argc; i++){ //Tentativa de leitura de cada um dos arquivos passados como parâmetro
        /*Leitura do Arquivo*/
        char * nomeDoArquivo = argv[i];
        FILE *file;
        file = fopen(nomeDoArquivo, "r"); 
        if (file == NULL) { //Verifica se o caminho existe
            break;
        } else {
            //<-----Ler campos do arquivo e adicionar em cada uma das variáveis aqui ----->
            while((fscanf(file,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &tempo, &alpha, &beta, &x0, &y0, &z0, &r0, &xl0, &yl0, &zl0, &vi, &xf, &yf, &zf, &rf, &dxf, &dyf, &zf, &vf)) != EOF){ //Enquanto não for o fim do arquivo
                // Tempo Alpha Beta X0 y0 z0 r0 xl0 yl0 zl0 |Vi| xf yf zf rf dxf dyf dzf |Vf|
                // 456.000000 104 89 -0.725655 2.910444 0.052357 3.000000 0.005108 -0.006719 -0.000104 0.008441 0.000000 0.000000 0.000000 0.000000 -0.001749 -0.005737 -0.000121 0.005999
                raio = EARTH_RADIUS + r0;
                w = sqrt(MI/(raio*raio*raio));
                calcularRendezvousBisseccao(x0,y0,z0,xl0,yl0,zl0,w);


            }
        }
    }
}