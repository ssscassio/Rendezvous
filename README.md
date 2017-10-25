# Rendezvous

Projeto para encontrar as variáveis físicas capazes de tornar o Rendezvous (Encontro com posição e velocidade relativa iguais a 0) de um veículo espacial com um detrito em órbita possível.

## Organização de Pastas


### 1. Arquivos de Entrada
Esta pasta contém exemplos dos arquivos de entrada que devem ser utilizados para a execução do código Rendezvous bem como a sua formatação correta

#### Debug Serial
Esta pasta contém um exemplo de arquivo para executar o teste em modo debug do código Serial do Rendezvous. Este arquivo deve conter apenas uma linha com a seguinte ordem de termos:

```
Tempo Alpha Beta x0 y0 z0 r0 xl0 yl0 zl0 |Vi| xf yf zf rf dxf dyf dzf |Vf| gama chi ve vex vey vez
```

Os dois outros arquivos **'v-0.004-0.005.dat'** e **'v-500-entradas.dat'** representam respectivamente um arquivo de entrada real da pesquisa do Rendezvous e um arquivo de entrada limitado a 500 linhas de conjuntos de parâmetros 

### 2. Arquivos de Saída
Esta pasta contém alguns exemplos dos arquivos gerados pelas execuções dos códigos em seus diferentes modos

#### 2.1. Debug Serial
O arquivo gerado pela execução serial no modo Debug contém informações de cada um dos coeficientes gerados a partir dos dados de entrada para calcular a equação que define o Rendezvous.

O arquivo gerado tem o formato de separação por virgula (CSV) e seu nome é "a-output-b.csv", em que 'a' é o indice do arquivo de parâmetro iniciando pelo valor 1 e 'b' é o indice da linha do arquivo sendo lido iniciando pelo valor 0.

#### 2.2. Serial
As linhas iniciais indicam os valores de entrada (posição e velocidade relativa) que geraram os resultados apresentados nas linhas seguintes.

#### 2.3. Paralelo
O arquivo gerado é diferente do resultante via código serial, porque os cálculos são executados assícronamente, gerando as linhas de resultado fora de ordem. Estes arquivos não apresentam resultados em ordem e nem um cabeçalho com os dados de entrada.

Recomenda-se ter em mãos o arquivo de entrada e fazer a mesclagem da linha do arquivo de entrada com o nome do arquivo de saída.(Indice da linha = 'b' no nome a-output-b.csv)

### 3. Planilha de Testes
Uma planilha Excel que pode ser utilizada para validar os resultados gerados pelos códigos aqui apresentados. O funcionamento da planilha baseia-se em preencher os dados de entrada do Rendezvous e o conjunto das três configurações físicas (Gama, Chi, Ve). A partir da todos os valores intermediários da equações serão calculados.

Esse planilha é útil para verificar se o código Serial está executando corretamente, comparando o resultado da execução em modo Debug com as informações calculadas na planilha.

### 4. Projeto Serial
Esta pasta contém o código desenvolvido para a execução do projeto de forma serial.Também existem informações detelhadas sobre como compilar e executar o projeto.

### 5. Projeto Paralelo (CUDA)
Esta pasta contém os códigos para a execução do projeto na arquitetura CUDA. Também existem informações detelhadas sobre como compilar e executar o projeto.

