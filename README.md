### _Integrantes do grupo_
***
_Carlos Henrique Goebel Teixeira Furtado_

_Lorenzo Baldissera Saldanha_

_Rodrigo Iasculski Da Conceição_
***
### _IPPD-Genetica_
Trabalho feito para disciplina de Introdução Processamento Paralelo E Distribuído na UFPEL.

Use para compilar:
```bash
  mpicc -fopenmp -O2 TrabalhoGeneticaIPPD.c -o TrabalhoGeneticaIPPD
```
E use para executar:
```bash
  ./TrabalhoGeneticaIPPD
```
***
### _Usando o Xivoco_
Usamos 8 nós para a execução.

Faça isso tudo no master:
```bash
  mpicc -fopenmp TrabalhoGeneticaIPPD.c -o TrabalhoGeneticaIPPD -lm
```
Após isso:
```bash
  scp ./TrabalhoGeneticaIPPD worker-1:~/
  scp ./TrabalhoGeneticaIPPD worker-2:~/
  scp ./TrabalhoGeneticaIPPD worker-3:~/
  scp ./TrabalhoGeneticaIPPD worker-4:~/
  scp ./TrabalhoGeneticaIPPD worker-5:~/
  scp ./TrabalhoGeneticaIPPD worker-6:~/
  scp ./TrabalhoGeneticaIPPD worker-7:~/
```
Logo após:
```bash
  mpirun -np 8 --host master,worker-1,worker-2,worker-3,worker-4,worker-5,worker-6,worker-7 ./TrabalhoGeneticaIPPD caenorhabditis_elegans.PRJNA13758.WBPS19.genomic.fa
```
