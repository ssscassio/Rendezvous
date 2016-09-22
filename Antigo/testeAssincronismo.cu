#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

/*
http://on-demand.gputechconf.com/gtc/2014/presentations/S4158-cuda-streams-best-practices-common-pitfalls.pdf

http://www.nvidia.com/content/gtc/documents/1122_gtc09.pdf


Possibilidade de utilizar CudaEvents para sinalizar a finalizar de uma serie de operações ou o espaço no buffer sendo utilizado.

**/
__global__ void somarVetores(double *A[], double*B[], double*C[]){

  int i = threadIdx.x;
  int j = threadIdx.y;

  C[i][j]=A[i][j]+B[i][j];
}

int main(){


  cudaStream_t stream1, stream2;
  cudaStreamCreate(&stream1);
  cudaStreamCreate(&stream2);

  double  A_h[10][10],B_h[10][10],C_h[10][10]; //Variaveis no host
  double  *A_d,*B_d,*C_d; //Variaveis do Device

  //Inicializando vetores
  for(int i = 0; i <10; i++){
    for(int j = 0; j <10; j++){
      A_h[i][j] = i*j;
      B_h[i][j] = i*j+1;
    }
  }

  dim3 threadsPerBlock(10,10);


  size_t size = 10*10*sizeof(double);

  cudaMalloc((void**)&A_d, size);
  cudaMalloc((void**)&B_d, size);
  cudaMalloc((void**)&C_d, size);

  cudaEvent_t start,stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  float gpu_time = 0.0f;

  cudaEventRecord(start, 0);
  cudaMemcpyAsync(&A_d, &A_h, size, cudaMemcpyHostToDevice,0);
  cudaMemcpyAsync(&B_d, &B_h, size, cudaMemcpyHostToDevice,0);
  cudaMemcpyAsync(&C_d, &C_h, size, cudaMemcpyHostToDevice,0);
  somarVetores<<<1,threadsPerBlock,0,0>>>(A_d,B_d,C_d);
  cudaMemcpyAsync(&C_h, &C_h, size, cudaMemcpyDeviceToHost,0);
  cudaEventRecord(stop, 0);

  while (cudaEventQuery(stop) == cudaErrorNotReady)
     {
         counter++;
     }

  // print the cpu and gpu times
  printf("time spent executing by the GPU: %.2f\n", gpu_time);
  printf("CPU executed %lu iterations while waiting for GPU to finish\n", counter);

      // release resources
  cudaEventDestroy(start);
  cudaEventDestroy(stop);
  cudaFreeHost(A_h);
  cudaFreeHost(B_h);
  cudaFreeHost(C_h);
  cudaFree(A_d);
  cudaFree(B_d);
  cudaFree(C_d);




}
