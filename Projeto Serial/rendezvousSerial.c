/**
    C++, rendezvousSerial.cu
    Purpose: Calcular as variáveis físicas de um veículo espacial capazes de tornar o Rendezvous possível

    @author Cássio Santos
    @version 1.0 31/07/17 
*/

/**
    Para compilaçao:
        ex.: gcc -o rendezvous rendezvousSerial.c -lm -std=c99
    Para execução:
        ex.: ./rendezvous v-0.005-0.006.dat v-0.007-0.008.dat
    Para debug:
        ex.: ./rendezvous debug debugInput.in > debugOutput.out

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//Declaração de constantes
#define MI 398600.4418
#define EARTH_RADIUS 6378.0
#define _USE_MATH_DEFINES
#define N_PRECISION 10

// Variável global para definir se a execução deve ser feita em modo de Debug ou não
// 0 - Falso - Debug Desativado
// 1 - Verdadeiro - Debug Ativado
int DEBUG = 0;

// Protótipos das funções para encontrar os coeficientes da equação do Rendezvous

double brute_A (int, double, double, double, double, double, double, double, double, double, int, double, double, double);

double brute_B (int, double, double, double, double, double, double, double, double, double, int, double, double, double);

double brute_C (int, double, double, double, double, double, double, double, double, double, int, double, double, double);

double brute_D (int, double, double, double, double, double, double, double, double, double, int, double, double, double);

double brute_E (int, double, double, double, double, double, double, double, double, double, int, double, double, double);

double brute_F (int, double, double, double, double, double, double, double, double, double, int, double, double, double);

double brute_G (int, double, double, double, double, double, double, double, double, double, int, double, double, double);

double brute_H (int, double, double, double, double, double, double, double, double, double, int, double, double, double);

double brute_I (int, double, double, double, double, double, double, double, double, double, int, double, double, double);

double brute_J (int, double, double, double, double, double, double, double, double, double, int, double, double, double);

// Protótipo da função que calcula o valor dos A sufixados na equação do Rendezvous
void brute_all(double *, int, double, double, double, double, double, double, double, int, double, double, double, double);



double calcularDiferenca (int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double gama, double chi, double w, double vex, double vey, double vez){
    double A[13];
    double a,b,c,d;
    double a1,a2,a3,a4,a5,a6;
    double result;

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

    result = pow(b-d,2)-4*(a-c)*(b*c-a*d);
    
    if(DEBUG){
        printf("a1: %lf\n", a1);
        printf("a2: %lf\n", a2);
        printf("a3: %lf\n", a3);
        printf("a4: %lf\n", a4);
        printf("a5: %lf\n", a5);
        printf("a6: %lf\n", a6);
        printf("a: %lf\n", a);
        printf("b: %lf\n", b);
        printf("c: %lf\n", c);
        printf("d: %lf\n", d);
        printf("Result: %lf\n", result);
    }
    return result;
}

void calcularRendezvousBisseccao(FILE * file, double x0, double y0, double z0, double xl0, double yl0, double zl0, double w) {
    double vex,vey,vez; //Variaveis Iteraveis //Componentes das Velocidades de exaustao
    double chi; //Variavel Iteravel
    double gama; //Variavel Iteravel

    //Variaveis de iteração
    int i;
    int iteracoes;

    //Variaveis para iteração por bissecção
    //X diz respeito ao valor de Ve sendo Iterado
    //Y diz respeito ao resultado da diferença
    double xInicial, xFinal, xMedio, yMedio, yInicial, yFinal;

    fprintf(file, "x0, y0, z0, xl0, yl0, zl0, w,\n");
    fprintf(file, "%lf, %lf, %lf, %lf, %lf, %lf,\n\n", x0,y0,z0, xl0, yl0, zl0);
    fprintf(file, "gama, chi, vex, y,\n");

    for (i = -14; i <= 2 ; i++) {//Iterando gama
        gama = pow(10,i);
        for(chi = 1; chi <= 100 ; chi++) {//Iterando Chi
            for(vex = 0.1; vex <= 10 ; vex+=0.1) {//Iterando nas componentes das velocidades de exaustão
                xInicial = vex;
                xFinal = vex+0.1;
                
                //Calculando diferença no ponto xInicial
                vex = vey = vez = xInicial;
                yInicial = calcularDiferenca(N_PRECISION, x0, y0, z0, xl0, yl0, zl0, gama, chi, w, vex, vey, vez);

                fprintf(file, "%.14lf, %lf, %lf, %lf,\n", gama, chi, xInicial, yInicial);
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
    double r0; //Altitude 
    
    //Variaveis inutilizadas (Por enquanto)
    double tempo, alpha, beta, vi, xf, yf, zf, rf, dxf, dyf, vf;
    //Variaveis usadas em modo Debug
    double gama, chi, ve, vex, vey, vez;

    //Verificação se algum nome de arquivo foi passado como parâmetro para a execução do código
    if(argc == 1){
        printf("Passe o nome dos arquivos de input como parâmetro\n");
        return 1;
    }

    //Verificação se o primeiro parâmetro passado na execução é a string 'debug'
    //O que definirá a variável global como 1, ativando assim o modo de Debug
    if(strcmp(argv[1],"debug") == 0){
        DEBUG = 1;
    }

    // Execução no modo de Debug o qual espera um arquivo com linhas formatadas assim como
    // a execução normal  com a adição das variáveis físicas ao final no formato:
    // Tempo Alpha Beta X0 y0 z0 r0 xl0 yl0 zl0 |Vi| xf yf zf rf dxf dyf dzf |Vf| Gama Chi Ve vex vey vez
    // 456.000000 104 89 -0.725655 2.910444 0.052357 3.000000 0.005108 -0.006719 -0.000104 0.008441 0.000000 0.000000 0.000000 0.000000 -0.001749 -0.005737 -0.000121 0.005999 100.0 10.0 5.0 5.0 5.0 5.0
    // Porém só executará a primeira linha do arquivo (Caso tenha mais de uma)
    if(DEBUG){
        char * nomeDoArquivo = argv[2];
        FILE *file;
        file = fopen(nomeDoArquivo, "r");
        //Verifica se o caminho existe
        if (file == NULL) { 
            printf("O arquivo ou caminho especificado não pôde ser encontrado\n");
            return 1;
        } else {
            // Leitura da primeira linha do arquivo
            fscanf(file,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &tempo, &alpha, &beta, &x0, &y0, &z0, &r0, &xl0, &yl0, &zl0, &vi, &xf, &yf, &zf, &rf, &dxf, &dyf, &zf, &vf , &gama, &chi, &ve, &vex, &vey, &vez);
            // Calculo da variável W baseado nas contantes e no raio lido apartir do arquivo
            raio = EARTH_RADIUS + r0;
            w = sqrt(MI/(raio*raio*raio));

            printf("===Dados de Entrada===\n");
            printf(" N: %d\n x0: %lf\n y0: %lf\n z0: %lf\n xl0: %lf\n yl0: %lf\n zl0: %lf\n gama: %lf\n chi: %lf\n r0: %lf\n ve: %lf\n vex: %lf\n vey: %lf\n vez: %lf\n",10, x0, y0, z0, xl0, yl0, zl0, gama, chi, r0, ve, vex, vey, vez);
            printf("======================\n");
            calcularDiferenca(N_PRECISION, x0, y0, z0, xl0, yl0, zl0, gama, chi, w, vex, vey, vez);
            return 0;

        }
    }

    //Tentativa de leitura de cada um dos arquivos passados como parâmetro
    for(int i = 1; i < argc; i++){ 
        /*Leitura do Arquivo*/
        char * nomeDoArquivo = argv[i];
        FILE *file;
        file = fopen(nomeDoArquivo, "r");
        int b = 0; //Variável para indice da linha do arquivo

        //Verifica se o caminho existe
        if (file == NULL) { 
            break;
        } else {
            //<-----Ler campos do arquivo e adicionar em cada uma das variáveis ----->
            // Exemplo de formatação de linha esperada no arquivo de entrada:
            // Tempo Alpha Beta X0 y0 z0 r0 xl0 yl0 zl0 |Vi| xf yf zf rf dxf dyf dzf |Vf|
            // 456.000000 104 89 -0.725655 2.910444 0.052357 3.000000 0.005108 -0.006719 -0.000104 0.008441 0.000000 0.000000 0.000000 0.000000 -0.001749 -0.005737 -0.000121 0.005999
            while((fscanf(file,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &tempo, &alpha, &beta, &x0, &y0, &z0, &r0, &xl0, &yl0, &zl0, &vi, &xf, &yf, &zf, &rf, &dxf, &dyf, &zf, &vf)) != EOF){ //Enquanto não for o fim do arquivo
                char nomeDoArquivoDeEscrita[256];
                sprintf( nomeDoArquivoDeEscrita, "%d-output-%d.csv", i , b);
                // Calculo da variável W baseado nas contantes e no raio lido apartir do arquivo
                raio = EARTH_RADIUS + r0;
                w = sqrt(MI/(raio*raio*raio));
                
                FILE *fileToWrite;
                fileToWrite = fopen(nomeDoArquivoDeEscrita, "w");
                
                calcularRendezvousBisseccao(fileToWrite,x0,y0,z0,xl0,yl0,zl0,w);
                b++;
            }
        }
    }
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
* @param vex Variável física da Velocidade de exaustão no eixo Y a ser calculado os coeficientes A sufixados
* @param vex Variável física da Velocidade de exaustão no eixo Z a ser calculado os coeficientes A sufixados
*/
void brute_all(double *A, int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, int X, double w, double vex, double vey, double vez){

    double a, B, D, e, G, H, I;    
    //Calculando valores de A, B, C, D, E, F, G, H, I e J
    a = brute_A(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 0, vex, vey,vez); 
    B = brute_B(N,x0,y0,z0,xl0,yl0,zl0,  Y, X, w, 0, vex, vey,vez);
    D = brute_D(N,x0,y0,z0,xl0,yl0,zl0,Y, X, w, 0, vex, vey,vez);
    e = brute_E(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 0, vex, vey,vez);
    G = brute_G(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 0, vex, vey,vez);
    H = brute_H(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 0, vex, vey,vez);
    I = brute_I(N,x0,y0,z0,xl0,yl0,zl0,  Y, X, w, 0, vex, vey,vez);

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

    //Debug
    if(DEBUG){
        for(int i = 1; i < 13 ; i++){
            printf("A[%d]: %lf\n",i, A[i]);
        }
    }
}


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
* @param vex Variável física da Velocidade de exaustão no eixo Y a ser calculado o valor de A
* @param vex Variável física da Velocidade de exaustão no eixo Z a ser calculado o valor de A
* @returns O coeficiênte A dado os valores iniciais e as variáveis físicas a serem testadas
*/
double brute_A (int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, double X, double w, int a, double vex, double vey, double vez) {
    double result = 0;
    int n;
    double aux;
    double sum = 0;

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

    if(DEBUG){
        printf("A: %lf\n", result);
    }
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
* @param vex Variável física da Velocidade de exaustão no eixo Y a ser calculado o valor de B
* @param vex Variável física da Velocidade de exaustão no eixo Z a ser calculado o valor de B
* @returns O coeficiênte B dado os valores iniciais e as variáveis físicas a serem testadas
*/
double brute_B (int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, double X, double w, int a, double vex, double vey, double vez){
    double result = 0;
    double sum = 0;
    int n;
    double aux;

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

    if(DEBUG){
        printf("B: %lf\n", result);
    }
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
* @param vex Variável física da Velocidade de exaustão no eixo Y a ser calculado o valor de C
* @param vex Variável física da Velocidade de exaustão no eixo Z a ser calculado o valor de C
* @returns O somatório dos coeficiêntes Cn dado os valores iniciais e as variáveis físicas a serem testadas
*/
double brute_C (int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, double X, double w, int a, double vex, double vey, double vez){
    double result = 0;
    int n;
    double aux;

    //Calculo do somatorio Cn
    for (n = 1; n <= N; n++) {
        aux = pow(n,a)*vex/(n*pow(X,n)*(1+pow((n*Y/w),2))) + n*Y*pow(n,a)*vey/(n*pow(X,n)*pow(w,2)*(1+pow((n*Y/w),2)));
        if (n%2 == 0) {
            aux = -aux;
        }

        result +=aux;   
    }


    if(DEBUG){
        printf("Cn*n^%d: %lf\n",a, result);
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
* @param vex Variável física da Velocidade de exaustão no eixo Y a ser calculado o valor de D
* @param vex Variável física da Velocidade de exaustão no eixo Z a ser calculado o valor de D
* @returns O coeficiênte D dado os valores iniciais e as variáveis físicas a serem testadas
*/
double brute_D (int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, double X, double w, int a, double vex, double vey, double vez) {
    double result = 0;

    result -= (2*vex* log((X+1)/X))/w;
    result += 4*y0 - 2*xl0/w;

    if(DEBUG){
        printf("D: %lf\n", result);
    }
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
* @param vex Variável física da Velocidade de exaustão no eixo Y a ser calculado o valor de E
* @param vex Variável física da Velocidade de exaustão no eixo Z a ser calculado o valor de E
* @returns O coeficiênte E dado os valores iniciais e as variáveis físicas a serem testadas
*/
double brute_E (int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, double X, double w, int a, double vex, double vey, double vez) {
    double result = 0;

    result -=  3*vex*log((X+1)/X);
    result +=  6*w*y0 - 3*xl0;

    if(DEBUG){
        printf("E: %lf\n", result);
    }
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
* @param vex Variável física da Velocidade de exaustão no eixo Y a ser calculado o valor de Fn
* @param vex Variável física da Velocidade de exaustão no eixo Z a ser calculado o valor de Fn
* @returns O somatório coeficiênte Fn dado os valores iniciais e as variáveis físicas a serem testadas
*/
double brute_F(int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, double X, double w, int a, double vex, double vey, double vez) {
    double result = 0;
    double sum = 0;
    int n;
    double aux;


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
 

    if(DEBUG){
        printf("Fn*n^%d: %lf\n",a, result);
    }
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
* @param vex Variável física da Velocidade de exaustão no eixo Y a ser calculado o valor de G
* @param vex Variável física da Velocidade de exaustão no eixo Z a ser calculado o valor de G
* @returns O coeficiênte G dado os valores iniciais e as variáveis físicas a serem testadas
*/
double brute_G (int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, double X, double w, int a, double vex, double vey, double vez) {
    double result = 0;
    double sum = 0;
    int n;
    double aux;

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

    if(DEBUG){
        printf("G: %lf\n", result);
    }
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
* @param vex Variável física da Velocidade de exaustão no eixo Y a ser calculado o valor de H
* @param vex Variável física da Velocidade de exaustão no eixo Z a ser calculado o valor de H
* @returns O coeficiênte H dado os valores iniciais e as variáveis físicas a serem testadas
*/
double brute_H (int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, double X, double w, int a, double vex, double vey, double vez) {
    double result = 0;
    double sum = 0;
    int n;
    double aux;
    
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

    if(DEBUG){
        printf("H: %lf\n", result);
    }
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
* @param vex Variável física da Velocidade de exaustão no eixo Y a ser calculado o valor de I
* @param vex Variável física da Velocidade de exaustão no eixo Z a ser calculado o valor de I
* @returns O coeficiênte I dado os valores iniciais e as variáveis físicas a serem testadas
*/
double brute_I (int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, double X, double w, int a, double vex, double vey, double vez) {
    double result = 0;
    double sum = 0;
    int n;
    double aux;
    
    result = zl0/w - (vez/w)*(log((X+1)/X));

    //Calculo do somatorio
    for (n = 1; n <= N; n++) {
        aux = ((vez)/(pow(n,2)*pow(X,n)*w))/(1+pow((n*Y)/w,2));
        if (n%2 == 0) {
            aux = -aux;
        }
        sum += aux;
    }

    result += sum;

    if(DEBUG){
        printf("I: %lf\n", result);
    }
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
* @param vex Variável física da Velocidade de exaustão no eixo Y a ser calculado o valor de Jn
* @param vex Variável física da Velocidade de exaustão no eixo Z a ser calculado o valor de Jn
* @returns O somatório coeficiênte Jn dado os valores iniciais e as variáveis físicas a serem testadas
*/
double brute_J(int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, double X, double w, int a, double vex, double vey, double vez){
    double result = 0;
    double sum = 0;
    int n;
    double aux;

    for (n = 1; n <= N; n++) {
        aux = vez/(n*pow(X,n)*w)/(1+pow((n*Y)/w,2));
        if (n%2 == 0) {
            aux = - aux;
        }
        aux *= pow(n,a);
        sum += aux;
    }
    
    result = sum;
    
    if(DEBUG){
        printf("Jn*n^%d: %12.10f\n",a, result);
    }
    return result;
}