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

    __syncthreads();
}

double __device__ calcularDiferenca (int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double gama, double chi, double w, double vex, double vey, double vez){
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

    return result;
}

void __global__ calcularRendezvousDevice(double d_x0,double d_y0,double d_z0,double d_xl0,double d_yl0,double d_zl0,double d_w){

  double gama = blockIdx.x; //Y(Gama) recebe o valor x atual do bloco;
  double chi = threadIdx.x; //X(Chi) recebe o valor y atual do thread;
  double ve = blockIdx.y; //ve||vex(Velocidade de exaustão) recebe o valor z atual do thread;

  double yInicial = calcularDiferenca(10, d_x0, d_y0, d_z0, d_xl0, d_yl0, d_zl0, gama, chi, d_w, ve, ve, ve);


}

void calcularRendezvous(double x0, double y0, double z0, double xl0, double yl0, double zl0, double w) {

    int tx = 16; //Número de iterações em Gama (de -14 a 2 com passo 1)
    int ty = 99; //Número de iterações em ve (de 0.1 a 10 com passo 0.1)
    int tz = 99; //Número de iterações em Chi (de 1 a 100 com passo 1)
    dim3 numBlocks(tx,ty);
    size_t size = sizeof(double);
    dim3 threadsPerBlock(tz);

    //Alocando a memória das variáveis no Device
    cudaMalloc((void**)&d_x0, size);
    cudaMalloc((void**)&d_y0, size);
    cudaMalloc((void**)&d_z0, size);
    cudaMalloc((void**)&d_xl0, size);
    cudaMalloc((void**)&d_yl0, size);
    cudaMalloc((void**)&d_zl0, size);
    cudaMalloc((void**)&d_w, size);

    //Copiando as variáveis do host para o Device
    cudaMemcpy(&d_x0, &x0, size, cudaMemcpyHostToDevice);
    cudaMemcpy(&d_y0, &y0, size, cudaMemcpyHostToDevice);
    cudaMemcpy(&d_z0, &z0, size, cudaMemcpyHostToDevice);
    cudaMemcpy(&d_xl0, &xl0, size, cudaMemcpyHostToDevice);
    cudaMemcpy(&d_yl0, &yl0, size, cudaMemcpyHostToDevice);
    cudaMemcpy(&d_zl0, &zl0, size, cudaMemcpyHostToDevice);
    cudaMemcpy(&d_w, &w, size, cudaMemcpyHostToDevice);

    calcularRendezvousDevice<<<numBlocks, threadsPerBlock>>>(d_x0,d_y0,d_z0,d_xl0,d_yl0,d_zl0,d_w);

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
    //Variaveis para modo Debug
    double gama, chi, ve, vex, vey, vez;

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
            while((fscanf(file,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &tempo, &alpha, &beta, &x0, &y0, &z0, &r0, &xl0, &yl0, &zl0, &vi, &xf, &yf, &zf, &rf, &dxf, &dyf, &zf, &vf)) != EOF){ //Enquanto não for o fim do arquivo
                // Tempo Alpha Beta X0 y0 z0 r0 xl0 yl0 zl0 |Vi| xf yf zf rf dxf dyf dzf |Vf|
                // 456.000000 104 89 -0.725655 2.910444 0.052357 3.000000 0.005108 -0.006719 -0.000104 0.008441 0.000000 0.000000 0.000000 0.000000 -0.001749 -0.005737 -0.000121 0.005999
                raio = EARTH_RADIUS + r0;
                w = sqrt(MI/(raio*raio*raio));
                calcularRendezvous(x0,y0,z0,xl0,yl0,zl0,w);

            }
        }
    }

}