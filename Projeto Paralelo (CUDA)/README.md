# Rendezvous CUDA

Projeto para Encontrar as variáveis físicas capazes de tornar o Rendezvous (Encontro com posição e velocidade relativa iguais a 0) de um veiculo espacial com um detrito em orbita possível.

## Getting Started

Estas instruções permitirão, caso tenha previamente o conjuto de variáveis que tornem o encontro possível, encontrar o conjunto de variáveis físicas que tornam o rendezvous possível para cada um desses conjutos de variáveis de entrada.

### Compilando

Para compilar o projeto:

```
nvcc rendezvousParalel.cu -o rendezvous
```

### Executando o projeto

Necessita de pelo menos um arquivo de entrada *input.dat* no qual cada uma das linhas representa um conjunto de variáveis de entrada das quais são respectivamente:
```
Tempo Alpha Beta X0 y0 z0 r0 xl0 yl0 zl0 |Vi| xf yf zf rf dxf dyf dzf |Vf|
``` 

Os arquivos de saida irão conter icognitas físicas descobertas bem como o resultado da equação.
 
Para executar o código execute o comando:
```
./rendezvous input1.dat input2.dat
```

### Resultados

Para cada arquivo de entrada na execução serão gerados tantos arquivos de saida quanto forem os números de linhas no arquivo

Por exemplo, para a execução do comando:
```
./rendezvous input1.dat input2.dat

```
Serão gerados arquivos no formato csv (Separados por virgula) com nomes:

```
"1-output-0.csv"
"1-output-0.csv"
.

"%d-output-%d.csv"
.
.
"2-output-500.csv"
```
Nos quais os indices representam respectivamente o arquivo de entrada e a linha do arquivo

Cada linha dos arquivos de saida estará no seguinte formato:

```
Gama , Chi , Ve , Diferenca
```

Na qual **Diferenca** representa o resultado da equação: *powf(b-d,2)-4*(a-c)*(b*c-a*d)*