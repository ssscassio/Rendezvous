#include <stdio.h>
#include <stdlib.h>
#include <math.h>

float funcao(double gama, double chi, double ve){
    return gama - chi + ve; //Expressão da equação aqui
}

int main(){

//Setar Gama, Setar Chi
//Definir intervalo inicial com base na variável restante(Ve)
//Definir Quantidade de intervalos que serão buscados (para ser paralelizado)
//
  int i, chi;
  float gama, vex;
  for (i = -14; i <=2; i++){
      gama = pow(10,i);
      for(chi = 1 ; chi<=100 ; chi++){
          float x0,x1,xm, valorxm;
          x0 = 0.1;
          x1 = 10;
          while(funcao(gama,chi,x0)*funcao(gama,chi,x1) <=0){
              xm = (x0+x1)/2; //Verificar isso pois pode sair do passo desejado
              valorxm = funcao(gama,chi,xm);
              if(abs(valorxm) < 0.0001){ //Calcular valor da função no ponto medio(Colocar precisão aqui)
                printf("gama= %f chi= %d ve= %f; funcao = %f\n",gama,chi,xm, funcao(gama,chi,xm));
                break;
              }
              if(funcao(gama,chi,xm)* funcao(gama,chi,x0) > 0){
                x0 =xm;
              }else{
                x1 =xm;
              }
          }
      }
    }
}
