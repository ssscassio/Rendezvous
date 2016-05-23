#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double ln(double x){
    double result;
    result = log(x)/log(M_E);
    return result;
}

//Método para achar de forma bruta o coeficiente A
double brute_A(int N, float x0, float y0, float z0, float xl0, float yl0, float zl0, float Y, int X, float w, int a, float vex, float vey, float vez){

double result =0;

//Calculo do somatorio
    int n;
     for(n = 1; n <=N; n++){
        if(n%2 ==0){
            result-= 2*vex/(n*pow(X,n)*w*(1+pow((n*Y/w),2)));
            result -= (n*Y*vey)/(n*pow(X,n)*pow(w,2)*(1+pow((n*Y/w),2)));
        }else{
            result+= 2*vex/(n*pow(X,n)*(1+pow((n*Y/w),2)));
            result += (n*Y*vey)/(n*pow(X,n)*pow(w,2)*(1+pow((n*Y/w),2)));
        }
    }

    result += 2*xl0/w - 3*y0;
    //Por mudança de base Ln((X+1)/X) = log((X+1)/X)/log(euler)
    result +=(2*vex*ln((X+1)/X))/w;

    return result;
}

// Método para achar de forma bruta o coeficiente B
double brute_B(int N, float x0, float y0, float z0, float xl0, float yl0, float zl0, float Y, int X, float w, int a, float vex, float vey, float vez){

double result =0;

//Calculo do somatorio
    int n;
    for(n = 1; n <=N; n++){
        if(n%2 ==0){
            result+= 2*vex*n*Y/(n*pow(X,n)*pow(w,2)*(1+pow((n*Y/w),2)));

            result += vey/(n*pow(X,n)*w*(1+pow((n*Y/w),2)));
        }else{
            result -= 2*n*vex*Y/(n*pow(X,n)*pow(w,2)*(1+pow((n*Y/w),2)));

            result -= vey/(n*pow(X,n)*w*(1+pow((n*Y/w),2)));
        }

    }


    result += yl0/w;
    //Por mudança de base Ln((X+1)/X) = log((X+1)/X)/log(euler)
    result += vey*(ln((X+1)/X))/w;


    return result;
}

// Método para achar de forma bruta o coeficiente C
double brute_C(int N, float x0, float y0, float z0, float xl0, float yl0, float zl0, float Y, int X, float w, int a, float vex, float vey, float vez){
/**
Usar variavel de entrada para definir se vai ou não multiplicar o N no somatório
pow(n,a)
*/
double result =0;

//Calculo do somatorio
   int n;
    for(n = 1; n <=N; n++){
        if(n%2 ==0){
            result += pow(n,a)*vex/(n*pow(X,n)*(1+pow((n*Y/w),2)));

            result += n*Y*pow(n,a)*vey/(n*pow(X,n)*pow(w,2)*(1+pow((n*Y/w),2)));
        }else{
            result-= pow(n,a)*vex/(n*pow(X,n)*(1+pow((n*Y/w),2)));

            result -= n*Y*pow(n,a)*vey/(n*pow(X,n)*pow(w,2)*(1+pow((n*Y/w),2)));
        }

    }

    return result;
}

// Método para achar de forma bruta o coeficiente D
double brute_D(int N, float x0, float y0, float z0, float xl0, float yl0, float zl0, float Y, int X, float w, int a, float vex, float vey, float vez){

double result =0;

    result -= (2*vex* ln((X+1)/X))/w;

    result += 4*y0 - 2*xl0/w;

    return result;
}

// Método para achar de forma bruta o coeficiente E
double brute_E(int N, float x0, float y0, float z0, float xl0, float yl0, float zl0, float Y, int X, float w, int a, float vex, float vey, float vez){

double result =0;

//Calculo do somatorio
    result -=  3*vex*ln((X+1)/X);

    result +=  6*w*y0 - 3*xl0;

    return result;
}

// Método para achar de forma bruta o coeficiente F
double brute_F(int N, float x0, float y0, float z0, float xl0, float yl0, float zl0, float Y, int X, float w, int a, float vex, float vey, float vez){
/**
Usar variavel de entrada para definir se vai ou não multiplicar o N no somatório
pow(n,a)
*/
double result =0;

//Calculo do somatorio
     int n;
    for(n = 1; n <=N; n++){
        if(n%2 ==0){
            result += pow(n,a)*vex*4/(n*pow(X,n)*n*Y*(1+pow((n*Y/w),2)));

            result += 2*vey*pow(n,a)/(n*pow(X,n)*w*(1+pow((n*Y/w),2)));

        }else{
            result -= pow(n,a)*vex*4/(n*pow(X,n)*n*Y*(1+pow((n*Y/w),2)));

            result -= 2*vey*pow(n,a)/(n*pow(X,n)*w*(1+pow((n*Y/w),2)));
        }

    }

    result -= vex/n*Y;

    return result;
}

// Método para achar de forma bruta o coeficiente G
double brute_G(int N, float x0, float y0, float z0, float xl0, float yl0, float zl0, float Y, int X, float w, int a, float vex, float vey, float vez){

double result =0;

//Calculo do somatorio
     int n;
    for(n = 1; n <=N; n++){
        if(n%2 ==0){
            result -= 3*vex/(pow(n,2)*pow(X,n)*w);
        }else{
            result += 3*vex/(pow(n,2)*pow(X,n)*w);
        }

    }


    result+= 2*yl0/w + x0;
    //Por mudança de base Ln((X+1)/X) = log((X+1)/X)/log(euler)
    result+= 2*vey*(ln((X+1)/X))/w;


    return result;
}

// Método para achar de forma bruta o coeficiente H
double brute_H(int N, float x0, float y0, float z0, float xl0, float yl0, float zl0, float Y, int X, float w, int a, float vex, float vey, float vez){

double result =0;

//Calculo do somatorio
     int n;
    for(n = 1; n <=N; n++){
        if(n%2 ==0){
            result += vez*Y/(pow(w,2)*pow(X,n)*(1+pow((n*Y/w),2)));
        }else{
            result -= vez*Y/(pow(w,2)*pow(X,n)*(1+pow((n*Y/w),2)));
        }

    }

    result+= z0;


    return result;
}

// Método para achar de forma bruta o coeficiente I
double brute_I(int N, float x0, float y0, float z0, float xl0, float yl0, float zl0, float Y, int X, float w, int a, float vex, float vey, float vez){

double result =0;

//Calculo do somatorio
     int n;
    for(n = 1; n <=N; n++){
        if(n%2 ==0){
            result += vez/(pow(n,2)*pow(X,n)*w*(1+pow((n*Y/w),2)));
        }else{
            result -= vez/(pow(n,2)*pow(X,n)*w*(1+pow((n*Y/w),2)));
        }
    }


    result+= zl0/w;
    //Por mudança de base Ln((X+1)/X) = log((X+1)/X)/log(euler)
    result -= vez*(ln((X+1)/X))/w;


    return result;
}

// Método para achar de forma bruta o coeficiente J
double brute_J(int N, float x0, float y0, float z0, float xl0, float yl0, float zl0, float Y, int X, float w, int a, float vex, float vey, float vez){
/**
Usar variavel de entrada para definir se vai ou não multiplicar o N no somatório
pow(n,a)
*/
double result =0;

//Calculo do somatorio
    int n;
    for(n = 1; n <=N; n++){
        if(n%2 ==0){
            result += vez*pow(n,a)/(n*pow(X,n)*(1+pow((n*Y/w),2)));

          }else{
            result -= vez*pow(n,a)/(n*pow(X,n)*(1+pow((n*Y/w),2)));
      }
    }

    return result;
}


double brute_all(double *A, int N, float x0, float y0, float z0, float xl0, float yl0, float zl0, float Y, int X, float w, float vex, float vey, float vez){


    double result = 0;

     //Calculando inicialmente a soma (A1² + A3² + A5²)
    A[1]= 2*brute_A(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 0, vex, vey,vez)*w + brute_E(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 0, vex, vey,vez) - Y*brute_F(N,x0,y0,z0, xl0,yl0, zl0, Y, X, w, 1, vex, vey,vez);

    A[3]= brute_B(N,x0,y0,z0,xl0,yl0,zl0,  Y, X, w, 0, vex, vey,vez)*w  - Y*brute_F(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 0, vex, vey,vez);
    A[5]= brute_I(N,x0,y0,z0,xl0,yl0,zl0,  Y, X, w, 0, vex, vey,vez)*w  - Y*brute_F(N,x0,y0,z0,xl0,yl0,zl0,  Y, X, w, 0, vex, vey,vez)*w+Y*brute_J(N,x0,y0,z0,xl0, yl0, zl0, Y, X, w, 0, vex, vey,vez)*w  - Y*brute_F(N,x0,y0,z0,xl0,yl0, zl0, Y, X, w, 1, vex, vey,vez);;

     //Calculando inicialmente a soma (A2² + A4² + A6²)
    A[2]= brute_G(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 0, vex, vey,vez) - 2*brute_B(N,x0,y0,z0,xl0,yl0,zl0,  Y, X, w, 0, vex, vey,vez) + brute_F(N,x0,y0,z0,xl0,yl0, zl0, Y, X, w, 0, vex, vey,vez);
    A[4]= brute_A(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 0, vex, vey,vez) +brute_D(N,x0,y0,z0,xl0,yl0,zl0,Y, X, w, 0, vex, vey,vez) + brute_C(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 0, vex, vey,vez);
    A[6]= brute_H(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 0, vex, vey,vez) - brute_J(N,x0,y0,z0,xl0,yl0,zl0,  Y, X, w, 0, vex, vey,vez);

     //Calculando inicialmente a soma (A7² + A9² + A11²)
    A[7]= 2*brute_B(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 0, vex, vey,vez)*w*w + Y*Y*brute_F(N,x0,y0,z0, xl0, yl0, zl0, Y, X, w, 2, vex, vey,vez);
    A[9]= Y*Y*brute_C(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 2, vex, vey,vez) - brute_A(N,x0,y0,z0,xl0,yl0,zl0,  Y, X, w, 0, vex, vey,vez)*w*w;
    A[11]= -w*w*brute_H(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 0, vex, vey,vez) - Y*Y*brute_J(N,x0,y0,z0,xl0,yl0,zl0, Y, X, w, 2, vex, vey,vez);

     //Calculando inicialmente a soma (A8² + A10² + A12²) A8=A1 // A10 = A3 // A12 = A5
    A[8]= A[1];
    A[10]=A[3];
    A[12]= A[5];


}

int main(){
    double A[13];
    float  gama,x0,y0,z0,xl0,yl0,zl0,vex,vey,vez,w,ve;
    printf("Secione o caso que desenha buscar:\n1) Encontrar ve para Y e X conhecidos\n2) Encontrar Y para ve e X conhecidos\n3) Encontrar X para Y e ve conhecidos\n");
    int selecao;
    FILE *arq;

    do{
      scanf("%d",&selecao);
    }while(selecao !=1 && selecao!=2 && selecao!=3);

    printf("Digite os valores de x0, y0, z0, xl0, yl0, zl0 e w respectivamente:\n");
    scanf("%f %f %f %f %f %f %f", &x0, &y0, &z0, &xl0, &yl0 ,&zl0, &w);
    float a,b,c,d;
    float a1,a2,a3,a4,a5,a6;
    int chi;
    int i;

    switch(selecao){

      case 1:
      //Caso 1

      arq = fopen("caso1.txt", "w");
      fprintf(arq, "Valores fixos: x0 = %f; y0 = %f; z0 = %f; xl0 = %f; yl0 = %f; zl0 = %f; w = %f\n\n\n",x0, y0, z0, xl0, yl0 ,zl0, w);
      for (i = -14; i <=2; i++){
          gama = pow(10,i);
          for(chi = 1 ; chi<=100 ; chi++){
              for(vex = 0.1; vex<=10; vex+=0.1){
                  vey =vez =vex;
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
                  float function = pow(b-d,2)-4*(a-c)*(b*c-a*d);
                  printf("Valor encontrado %f\n",function);
                  if( function >= -1 && function <=1){
                      ve = pow(3,1/2)*vex;
                      fprintf(arq, "Gama: 10^%d     Chi: %d       ve: %f        diff: %f\n", i,chi,ve, function);
                  }

              }}}

      break;
      case 2:
      //Caso 2
      arq = fopen("caso2.txt", "w");
      fprintf(arq, "Valores fixos: x0 = %f; y0 = %f; z0 = %f; xl0 = %f; yl0 = %f; zl0 = %f; w = %f\n\n\n",x0, y0, z0, xl0, yl0 ,zl0, w);
      for (ve = 0.1; ve <=10; ve+=0.1){
          for(chi = 1 ; chi<=100 ; chi++){
              for(i = -14; i<=2; i++){
                  vey =vez =vex = ve/pow(3,1/2);
                  gama = pow(10,i);
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
                  float function = pow(b-d,2)-4*(a-c)*(b*c-a*d);
                  printf("Valor encontrado %f\n",function);
                  if( function >= -1 && function <=1){
                    fprintf(arq, "ve: %f     Chi: %d       Gama: 10^%d        diff: %f\n", ve,chi, i, function);
                  }
              }
            }
          }

        break;
        case 3:
        //Caso 3
        arq = fopen("caso2.txt", "w");
      fprintf(arq, "Valores fixos: x0 = %f; y0 = %f; z0 = %f; xl0 = %f; yl0 = %f; zl0 = %f; w = %f\n\n\n",x0, y0, z0, xl0, yl0 ,zl0, w);
        for (ve = 0.1; ve <=10; i++){
            gama = pow(10,i);
              for(vex = 0.1; vex<10; vex+=0.1){
                for(chi = 1 ; chi<=100 ; chi++){
                    vey = vez =vex = ve/pow(3,1/2);
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
                    float function = pow(b-d,2)-4*(a-c)*(b*c-a*d);
                    printf("Valor encontrado %f\n",function);
                    if( function >= -1 && function <=1){

                        fprintf(arq, "Gama: 10^%d     Chi: %d       ve: %f        diff: %f\n", i,chi,ve, function);

                    }
                }
              }
            }

          break;

    }

    fclose(arq);

    /*while(1){
      float  gama,x0,y0,z0,xl0,yl0,zl0,vex,vey,vez,w,Y,X;
      printf("Digite os valores de x0, y0, z0, xl0,yl0, zl0, w,Y e X respectivamente:\n");
      scanf("%f %f %f %f %f %f %f %f %f", &x0, &y0, &z0, &xl0, &yl0 ,&zl0, &w ,&Y, &X);

    double a = brute_A(10, x0, y0, z0, xl0, yl0, zl0, Y, X,  w, a, vex, vey,  vez);
    printf("%.17lf\n",a );
  }*/
}
