# Rendezvous

Projeto para encontrar as variáveis físicas capazes de tornar o Rendezvous (Encontro com posição e velocidade relativa iguais a 0) de um veiculo espacial com um detrito em orbita possível.

## Organização de Pastas


### 1. Arquivos de Entrada
Essa pasta contém exemplos dos arquivos de entrada que devem ser utilizados para a execução do código Rendezvous bem como a sua formatação correta

#### Debug Serial
    Esta pasta contém um exemplo de arquivo para executar o teste em modo debug do código Serial do Rendezvous. Este arquivo deve conter apenas uma linha com a seguinte ordem de termos:

```
Tempo Alpha Beta X0 y0 z0 r0 xl0 yl0 zl0 |Vi| xf yf zf rf dxf dyf dzf |Vf| gama chi ve vex vey vez
```
Os dois outros arquivos **'v-0.004-0.005.dat'** e **'v-500-entradas.dat'** representam respectivamente um arquivo de entrada real da pesquisa do Rendezvous e um arquivo de entrada limitado a 500 linhas de conjuntos de parâmetros 

### 2. Arquivos de Saída
Esta pasta contém alguns exemplos dos arquivos gerados pelas execuções dos códigos em seus diferentes modos
#### Debug Serial
O arquivo gerado pela execução serial no modo Debug contém informações de cada um dos coeficientes gerados a partir dos dados de entrada para realizar a equação que define o Rendezvous
#### Serial
O nome do arquivo de saída é no formato a-output-b.csv onde a é o indice do arquivo de parâmetro iniciando pelo valor 1 e b é o indice da linha do arquivo sendo lido iniciando pelo valor 0.

O arquivo gerado tem o formato de separação por virgula (Comma Separated .csv) e apresenta nas suas linhas iniciais um cabeçalho indicando os valores de entrada (posição e velocidade relativa) que tal este resultado

#### Paralelo
O nome do arquivo de saída é no formato a-output-b.csv onde a é o indice do arquivo de parâmetro iniciando pelo valor 1 e b é o indice da linha do arquivo sendo lido iniciando pelo valor 0.

O arquivo gerado tem o formato de separação por virgula (Comma Separated .csv) e diferente do código serial, devido a execução dos calculos serem feitas de forma paralela e portanto executando os prints de forma randomica, os arquivos gerados pela execução paralela não apresentam uma ordem e nem um cabeçalho com os dados de entrada para essa execução

Recomenda-se ter em mãos o arquivo de entrada e fazer a mesclagem da linha do arquivo de entrada com o nome do arquivo de saída.(Indice da linha = b em a-output-b.csv)


### 3. Planilha de Testes
Essa pasta contém um arquivo de planilha no formato .xsls e pode ser utilizado para validar os resultados gerados pelos códigos aqui apresentados. A sua estrutra baseia-se em preencher os dados do arquivo de entrada do Rendezvous e escolher para qual conjunto das três configurações físicas (Gama, Chi, Ve) você deseja verificar passo a passo o desenvolvimento da equação do Rendezvous.

Você pode verificar se tudo está correto no código Serial por exemplo, comparando o resultado da execução em modo Debug com as informações apresentadas nessa planilha

### 4. Projeto Paralelo (CUDA)
Esta pasta contém os códigos desenvolvidos e necessários para a execução do projeto na arquitetura CUDA.
Nessa pasta podem ser encontradas informações mais detelhadas de como compilar e executar o projeto
### 5. Projeto Serial
Esta pasta contém o códigod desenvolvido para a execução do projeto de forma serial.
Nessa pasta podem ser encontradas informações mais detelhadas de como compilar e executar o projeto