# Rendezvous Serial

Projeto para Encontrar as variáveis físicas capazes de tornar o Rendezvous (Encontro com posição e velocidade relativa iguais a 0) de um veiculo espacial com um detrito em orbita possível.

## Getting Started

Estas instruções permitirão, caso tenha previamente o conjuto de variáveis que tornem o encontro possível, encontrar o conjunto de variáveis físicas que tornam o rendezvous possível para cada um desses conjutos de variáveis de entrada.

### Compilando

Para compilar o projeto:

```
gcc -o rendezvous RendezvousZeroFinding.c -lm -std=c99 -Wall -O3
```

### Executando o projeto

Necessita de pelo menos um arquivo de entrada *input.dat* no qual cada uma das linhas representa um conjunto de variáveis de entrada das quais são respectivamente:
```
Tempo Alpha Beta X0 y0 z0 r0 xl0 yl0 zl0 |Vi| xf yf zf rf dxf dyf dzf |Vf|
``` 

O arquivo de saida irá conter uma linha inicial contendo os valores de entrada e as linhas subsequentes contendo os conjuntos de icognitas físicas descobertas que satisfazem a equação.
 
Para executar o código execute o comando:
```
./rendezvous input1.dat input2.dat >output.dat

```

### Modo Debug
Para executar o código em modo Debug e encontrar o valor da função de Rendezvous para apenas 1 conjunto de dados especifivo, necessita de um arquivo com os seguintes conjuntos de dados.
Dispostos em uma unica linha do arquivo:
```
Tempo Alpha Beta X0 y0 z0 r0 xl0 yl0 zl0 |Vi| xf yf zf rf dxf dyf dzf |Vf| gama chi ve vex vey vez
```

o arquivo de saida conterá linhas correspondentes a cada uma das variaveis durante o processo de Debug as quais podem ser comparadas com os resultados obtidos através do arquivo `.xlsx` dentro do caminho `Projeto/Teste`

Para executar o código em modo debug execute o comando:
```
./rendezvous debug debugInput.in > debug.out
```
Sendo o parâmetro *debug* obrigatório como segundo parâmetro da execução e o parâmetro seguinte sendo o arquivo caminho do arquivo de dados de entrada.
