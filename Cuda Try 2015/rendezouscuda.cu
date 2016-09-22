#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

__device__ double ln(double x){
    return log(x)/log(M_E);
}


//Método para achar de forma bruta o coeficiente A
double __device__ brute_A(int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, int X, double w, int a, double vex, double vey, double vez){

  double result = 0;

  //Calculo do somatorio
  int n;
  for(n = 1; n <=N; n++){
    if(n%2==0){
      result -= (2*vex)/(n*powf(X,n)*w*(1+((n*Y)/w)*((n*Y)/w)));
      result -= (n*Y*vey)/(n*powf(X,n)*(w*w)*(1+(n*Y/w)*(n*Y/w)));

    }else{
      result += (2*vex)/(n*powf(X,n)*(1+(n*Y/w)*(n*Y/w)));
      result += (n*Y*vey)/(n*powf(X,n)*(w*w)*(1+(n*Y/w)*(n*Y/w)));
    }
  }

  result += (2*xl0)/w - 3*y0;
  //Por mudança de base Ln((X+1)/X) = log((X+1)/X)/log(euler)
  result +=(2*vex*ln((X+1)/X))/w;

  return result;
}

// Método para achar de forma bruta o coeficiente B
double __device__ brute_B(int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, int X, double w, int a, double vex, double vey, double vez){

  double result =0;

//Calculo do somatorio
  int n;
  for(n = 1; n <=N; n++){
    if(n%2 ==0){
        result+= (2*vex*n*Y)/(n*powf(X,n)*(w*w)*(1+((n*Y/w)*(n*Y/w))));

        result += vey/(n*powf(X,n)*w*(1+((n*Y/w)*(n*Y/w))));
    }else{
        result -= (2*n*vex*Y)/(n*powf(X,n)*powf(w,2)*(1+((n*Y/w)*(n*Y/w))));

        result -= vey/(n*powf(X,n)*w*(1+((n*Y/w)*(n*Y/w))));
    }
  }

  result += yl0/w;
  //Por mudança de base Ln((X+1)/X) = log((X+1)/X)/log(euler)
  result += (vey*(ln((X+1)/X)))/w;
  return result;
}

// Método para achar de forma bruta o coeficiente C
double __device__ brute_C(int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, int X, double w, int a, double vex, double vey, double vez){
  /**
  Usar variavel de entrada para definir se vai ou não multiplicar o N no somatório
  powf(n,a)
  */
  double result =0;

  //Calculo do somatorio
  int n;
  for(n = 1; n <=N; n++){
      if(n%2 ==0){
          result += (powf(n,a)*vex)/(n*powf(X,n)*(1+(n*Y/w)*(n*Y/w)));

          result += (n*Y*powf(n,a)*vey)/(n*powf(X,n)*(w*w)*(1+(n*Y/w)*(n*Y/w)));
      }else{
          result-= powf(n,a)*vex/(n*powf(X,n)*(1+powf((n*Y/w),2)));

          result -= n*Y*powf(n,a)*vey/(n*powf(X,n)*(w*w)*(1+powf((n*Y/w),2)));
      }

  }

  return result;
}

// Método para achar de forma bruta o coeficiente D
double __device__ brute_D(int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, int X, double w, int a, double vex, double vey, double vez){
  return ((4*y0) - ((2*xl0)/w) - ((2*vex*ln((X+1)/X))/w) );
}

// Método para achar de forma bruta o coeficiente E
double brute_E(int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, int X, double w, int a, double vex, double vey, double vez){
  return  6*w*y0 - 3*xl0 - 3*vex*ln((X+1)/X);
}

// Método para achar de forma bruta o coeficiente F
double __device__ brute_F(int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, int X, double w, int a, double vex, double vey, double vez){
  /**
  Usar variavel de entrada para definir se vai ou não multiplicar o N no somatório
  powf(n,a)
  */
  double result =0;

  //Calculo do somatorio
  int n;
  for(n = 1; n <=N; n++){
    if(n%2 ==0){
        result += (powf(n,a)*vex*4)/(n*powf(X,n)*n*Y*(1+((n*Y/w))*((n*Y)/w)));
        result += (2*vey*powf(n,a))/(n*powf(X,n)*w*(1+(n*Y/w)*(n*Y/w)));

    }else{
        result -= (powf(n,a)*vex*4)/(n*powf(X,n)*n*Y*(1+(n*Y/w)*(n*Y/w)));
        result -= (2*vey*powf(n,a))/(n*powf(X,n)*w*(1+(n*Y/w)*(n*Y/w)));
    }
  }

  result -= vex/(n*Y);
  return result;
}

// Método para achar de forma bruta o coeficiente G
double __device__ brute_G(int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, int X, double w, int a, double vex, double vey, double vez){

  double result =0;

  //Calculo do somatorio
  int n;
  for(n = 1; n <=N; n++){
    if(n%2 ==0){
        result -= (3*vex)/((n*n)*powf(X,n)*w);
    }else{
        result += (3*vex)/((n*n)*powf(X,n)*w);
    }
  }

  result+= (2*yl0)/w + x0;
  //Por mudança de base Ln((X+1)/X) = log((X+1)/X)/log(euler)
  result+= (2*vey*(ln((X+1)/X)))/w;

  return result;

}

// Método para achar de forma bruta o coeficiente H
double __device__ brute_H(int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, int X, double w, int a, double vex, double vey, double vez){

  double result =0;
  //Calculo do somatorio
  int n;
  for(n = 1; n <=N; n++){
    if(n%2 ==0){
        result += (vez*Y)/((w*w)*powf(X,n)*(1+((n*Y/w)*(n*Y/w))));
    }else{
        result -= (vez*Y)/((w*w)*powf(X,n)*(1+((n*Y/w)*(n*Y/w))));
    }
  }

  result+= z0;
  return result;
}

// Método para achar de forma bruta o coeficiente I
double __device__ brute_I(int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, int X, double w, int a, double vex, double vey, double vez){

  double result =0;
  //Calculo do somatorio
  int n;
  for(n = 1; n <=N; n++){
    if(n%2 ==0){
        result += vez/((n*n)*powf(X,n)*w*(1+((n*Y/w)*(n*Y/w))));
    }else{
        result -= vez/((n*n)*powf(X,n)*w*(1+((n*Y/w)*(n*Y/w))));
    }
  }

  result+= zl0/w;
  //Por mudança de base Ln((X+1)/X) = log((X+1)/X)/log(euler)
  result -= (vez*(ln((X+1)/X)))/w;

  return result;
}

// Método para achar de forma bruta o coeficiente J
double __device__ brute_J(int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, int X, double w, int a, double vex, double vey, double vez){
  /**
  Usar variavel de entrada para definir se vai ou não multiplicar o N no somatório
  powf(n,a)
  */
  double result =0;

  //Calculo do somatorio
  int n;
  for(n = 1; n <=N; n++){
    if(n%2 ==0){
      result += (vez*powf(n,a))/(n*w*powf(X,n)*(1+((n*Y/w)*(n*Y/w))));
    }else{
      result -= (vez*powf(n,a))/(n*w*powf(X,n)*(1+((n*Y/w)*(n*Y/w))));
    }
  }

  return result;
}


void __device__ brute_all(double *A, int N, double x0, double y0, double z0, double xl0, double yl0, double zl0, double Y, int X, double w, double vex, double vey, double vez){


   //Calculando inicialmente a soma (A1² + A3² + A5²)
  A[1]= 2*brute_A(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 0, vex, vey,vez)*w + brute_E(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 0, vex, vey,vez) - Y*brute_F(N,x0,y0,z0, xl0,yl0, zl0, Y, X, w, 1, vex, vey,vez);

  A[3]= brute_B(N,x0,y0,z0,xl0,yl0,zl0,  Y, X, w, 0, vex, vey,vez)*w - Y*brute_C(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 1, vex, vey,vez);
  //Errado: Brute nC em vez de Brute nF - Verificar com o professor

  A[5]= brute_I(N,x0,y0,z0,xl0,yl0,zl0,  Y, X, w, 0, vex, vey,vez)*w  + Y*brute_J(N,x0,y0,z0,xl0, yl0, zl0, Y, X, w, 1, vex, vey,vez);

   //Calculando inicialmente a soma (A2² + A4² + A6²)
  A[2]= brute_G(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 0, vex, vey,vez) - 2*brute_B(N,x0,y0,z0,xl0,yl0,zl0,  Y, X, w, 0, vex, vey,vez) + brute_F(N,x0,y0,z0,xl0,yl0, zl0, Y, X, w, 0, vex, vey,vez);

  A[4]= brute_A(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 0, vex, vey,vez) + brute_D(N,x0,y0,z0,xl0,yl0,zl0,Y, X, w, 0, vex, vey,vez) + brute_C(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 0, vex, vey,vez);

  A[6]= brute_H(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 0, vex, vey,vez) - brute_J(N,x0,y0,z0,xl0,yl0,zl0,  Y, X, w, 0, vex, vey,vez);

   //Calculando inicialmente a soma (A7² + A9² + A11²)
  A[7]= 2*brute_B(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 0, vex, vey,vez)*w*w + Y*Y*brute_F(N,x0,y0,z0, xl0, yl0, zl0, Y, X, w, 2, vex, vey,vez);

  A[9]= Y*Y*brute_C(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 2, vex, vey,vez) - brute_A(N,x0,y0,z0,xl0,yl0,zl0,  Y, X, w, 0, vex, vey,vez)*w*w;

  A[11]= -w*w*brute_H(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 0, vex, vey,vez) - Y*Y*brute_J(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 2, vex, vey,vez);

   //Calculando inicialmente a soma (A8² + A10² + A12²) A8=A1 // A10 = A3 // A12 = A5
  A[8]= A[1];

  A[10]=A[3];

  A[12]= A[5];

  __syncthreads();


}


__global__ void case1(double x0, double y0, double z0, double xl0, double yl0, double zl0, double w){

  FILE *arq;
  arq = fopen("caso1.txt", "w");


  double A[13];
  int N = 10;
  double a,b,c,d;
  double a1,a2,a3,a4,a5,a6;

  double Y = blockIdx.x; //Y(Gama) recebe o valor x atual do bloco;
  double X = threadIdx.x; //X(Chi) recebe o valor y atual do thread;
  double vex = threadIdx.y; //ve||vex(Velocidade de exaustão) recebe o valor z atual do thread;
  double gama, vey, vez,ve;

  //Tranformação dos Id atuais dos threads nos valores corretos.
  Y -=14;
  X += 1;
  vex += 1;
  vex *= 0.1;
  gama = powf(10,Y);

  //Calculando os valores dde Asub(n) e os valores de a b c e d
  vey = vez =vex;
  brute_all(A, N, x0, y0, z0, xl0, yl0, zl0, gama, X, w, vex, vey, vez);
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
  double function = (b-d)*(b-d)-4*(a-c)*(b*c-a*d);
  if( function >= -0.01 && function <=0.01){
      ve = powf(3,1/2)*vex;
      fprintf(arq, "Gama: 10^%d     Chi: %d       ve: %f        diff: %f\n", Y,X,ve, function);
  }

}

__global__ void case2(double x0, double y0, double z0, double xl0, double yl0, double zl0, double w){

  FILE *arq;
  arq = fopen("caso2.txt", "w");


  double A[13];
  int N = 10;
  double a,b,c,d;
  double a1,a2,a3,a4,a5,a6;
  double gama, vex, vey, vez;

  double Y = blockIdx.x; //Y(Gama) recebe o valor x atual do bloco;
  double X = threadIdx.x; //X(Chi) recebe o valor y atual do thread;
  double ve = threadIdx.y; //ve||vex(Velocidade de exaustão) recebe o valor z atual do thread;

  //Tranformação dos Id atuais dos threads nos valores corretos.
  Y -=14;
  X += 1;
  ve += 1;
  ve *= 0.1;

  //Fazendo devidas transformações
  gama = powf(10,Y);
  vey =vez =vex = ve/powf(3,1/2);

  //Calculando os valores dde Asub(n) e os valores de a b c e d

  brute_all(A, N, x0, y0, z0, xl0, yl0, zl0, gama, X, w, vex, vey, vez);
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
  double function = (b-d)*(b-d)-4*(a-c)*(b*c-a*d);
  if( function >= -0.01 && function <=0.01){
    fprintf(arq, "ve: %f     Chi: %d       Gama: 10^%d        diff: %f\n", ve,X, Y, function);
  }

}

__global__ void case3(double x0, double y0, double z0, double xl0, double yl0, double zl0, double w){

  FILE *arq;
  arq = fopen("caso3.txt", "w");


  double A[13];
  int N = 10;
  double a,b,c,d;
  double a1,a2,a3,a4,a5,a6;
  double gama, vex, vey, vez;

  double Y = blockIdx.x; //Y(Gama) recebe o valor x atual do bloco;
  double X = threadIdx.x; //X(Chi) recebe o valor y atual do thread;
  double ve = threadIdx.y; //ve||vex(Velocidade de exaustão) recebe o valor z atual do thread;

  //Tranformação dos Id atuais dos threads nos valores corretos.
  Y -=14;
  X += 1;
  ve += 1;
  ve *= 0.1;
  gama = powf(10,Y);
  vey = vez =vex = ve/powf(3,1/2);

  //Calculando os valores dde Asub(n) e os valores de a b c e d
  brute_all(A, N, x0, y0, z0, xl0, yl0, zl0, gama, X, w, vex, vey, vez);
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
  double function = (b-d)*(b-d)-4*(a-c)*(b*c-a*d);
  if( function >= -0.01 && function <=0.01){
      fprintf(arq, "Gama: 10^%d     ve: %f       Chi: %d        diff: %f\n", Y,ve,X, function);
  }

}

int main(){
  int selecao;
  double  x0,y0,z0,xl0,yl0,zl0,w;
  double  d_x0,d_y0,d_z0,d_xl0,d_yl0,d_zl0,d_w;

  FILE *arq;


  printf("Secione o caso que deseja buscar:\n1) Encontrar ve para Y e X conhecidos\n2) Encontrar Y para ve e X conhecidos\n3) Encontrar X para Y e ve conhecidos\n");
  do{
    scanf("%d",&selecao);
  }while(selecao !=1 && selecao!=2 && selecao!=3);

  printf("Digite os valores de x0, y0, z0, xl0, yl0, zl0 e w respectivamente:\n");
  scanf("%f %f %f %f %f %f %f", &x0, &y0, &z0, &xl0, &yl0 ,&zl0, &w);

  //int tx = 16; //Número de iterações em Gama (de -14 a 2 com passo 1)
  int ty = 99; //Número de iterações em Chi (de 1 a 100 com passo 1)
  int tz = 99; //Número de iterações em ve (de 0.1 a 10 com passo 0.1)
  dim3 numBlocks(16);
  size_t size = sizeof(double);
  dim3 threadsPerBlock(ty,tz);

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

//Colocar o switch aqui futuramente---
  switch(selecao){
  //Case 1:
  case 1:
    arq = fopen("caso1.txt", "w");
    fprintf(arq, "Valores fixos: x0 = %f; y0 = %f; z0 = %f; xl0 = %f; yl0 = %f; zl0 = %f; w = %f\n\n\n",x0, y0, z0, xl0, yl0 ,zl0, w);
    case1<<<numBlocks, threadsPerBlock>>>(d_x0,d_y0,d_z0,d_xl0,d_yl0,d_zl0,d_w);
    break;

  case 2:
    arq = fopen("caso2.txt", "w");
    fprintf(arq, "Valores fixos: x0 = %f; y0 = %f; z0 = %f; xl0 = %f; yl0 = %f; zl0 = %f; w = %f\n\n\n",x0, y0, z0, xl0, yl0 ,zl0, w);
    case2<<<numBlocks, threadsPerBlock>>>(d_x0,d_y0,d_z0,d_xl0,d_yl0,d_zl0,d_w);
    break;

  case 3:
    arq = fopen("caso3.txt", "w");
    fprintf(arq, "Valores fixos: x0 = %f; y0 = %f; z0 = %f; xl0 = %f; yl0 = %f; zl0 = %f; w = %f\n\n\n",x0, y0, z0, xl0, yl0 ,zl0, w);
    case3<<<numBlocks, threadsPerBlock>>>(d_x0,d_y0,d_z0,d_xl0,d_yl0,d_zl0,d_w);
    break;
  }
  //Liberando memória
  cudaFree(&d_x0); cudaFree(&d_y0); cudaFree(&d_z0); cudaFree(&d_xl0); cudaFree(&d_yl0); cudaFree(&d_zl0); cudaFree(&d_w);


}
