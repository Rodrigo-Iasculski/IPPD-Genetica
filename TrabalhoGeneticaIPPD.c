#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <omp.h>
#include <stdlib.h>

int main(int argc, char **argv) {
    int rank, size, provided;
    int tamanhoTotal = 0;
    int tamanhoSequencia, linhaTemp, linhaSequencia;
    int sequenciasMPI, resto, totalLocal;
    int offset, totaldeSend, totaldeRecv, idx, len;
    int ALocal = 0, CLocal = 0, TLocal = 0, GLocal = 0, GCLocal = 0;
    int ATotal = 0, CTotal = 0, TTotal = 0, GTotal = 0, GCTotal = 0;
    int totalBases = 0;
    double GCLocalAntes;
    FILE *arquivo = NULL;
    char *linha = NULL;
    char (*headersTotal)[512] = NULL;
    char **seqTotal = NULL;
    int *tamTotal = NULL;
    int *GCSTotal = NULL;
    int *capacidades = NULL;
    char **seqLocal = NULL;
    int *tamLocal = NULL;
    int *GCSLocal = NULL;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        arquivo = fopen("caenorhabditis_elegans.PRJNA13758.WBPS19.genomic.fa", "r");
        if (!arquivo) {
            printf("Erro ao abrir o arquivo\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        tamanhoSequencia = 10;
        linhaTemp = 65536;
        linhaSequencia = -1;

        linha = malloc(linhaTemp * sizeof(char));
        headersTotal = malloc(tamanhoSequencia * sizeof(*headersTotal));
        seqTotal = malloc(tamanhoSequencia * sizeof(char*));
        tamTotal = malloc(tamanhoSequencia * sizeof(int));
        GCSTotal = malloc(tamanhoSequencia * sizeof(int));
        capacidades = malloc(tamanhoSequencia * sizeof(int));

        while (fgets(linha, linhaTemp, arquivo)) {
            linha[strcspn(linha, "\r\n")] = 0;

            if (linha[0] == '>') {
                linhaSequencia++;
                if (linhaSequencia >= tamanhoSequencia) {
                    tamanhoSequencia *= 2;
                    headersTotal = realloc(headersTotal, tamanhoSequencia * sizeof(*headersTotal));
                    seqTotal = realloc(seqTotal, tamanhoSequencia * sizeof(char*));
                    tamTotal = realloc(tamTotal, tamanhoSequencia * sizeof(int));
                    GCSTotal = realloc(GCSTotal, tamanhoSequencia * sizeof(int));
                    capacidades = realloc(capacidades, tamanhoSequencia * sizeof(int));
                }

                strncpy(headersTotal[linhaSequencia], linha, 512 - 1);
                headersTotal[linhaSequencia][512 - 1] = '\0';

                capacidades[linhaSequencia] = 4096;
                seqTotal[linhaSequencia] = malloc(capacidades[linhaSequencia] * sizeof(char));
                seqTotal[linhaSequencia][0] = '\0';
                tamTotal[linhaSequencia] = 0;
                GCSTotal[linhaSequencia] = 0;

            } else if (linhaSequencia >= 0 && strlen(linha) > 0) {
                len = strlen(linha);

                while (tamTotal[linhaSequencia] + len + 1 > capacidades[linhaSequencia]) {
                    capacidades[linhaSequencia] *= 2;
                    seqTotal[linhaSequencia] = realloc(seqTotal[linhaSequencia], capacidades[linhaSequencia] * sizeof(char));
                }
                strcpy(seqTotal[linhaSequencia] + tamTotal[linhaSequencia], linha);
                tamTotal[linhaSequencia] += len;
            }
        }

        fclose(arquivo);
        free(linha);
        free(capacidades);
        tamanhoTotal = linhaSequencia + 1;
    }

    MPI_Bcast(&tamanhoTotal, 1, MPI_INT, 0, MPI_COMM_WORLD);

    sequenciasMPI = tamanhoTotal / size;
    resto = tamanhoTotal % size;
    if (rank == size - 1) {
        totalLocal = sequenciasMPI + resto;
    } else {
        totalLocal = sequenciasMPI;
    }

    seqLocal = malloc(totalLocal * sizeof(char*));
    tamLocal = malloc(totalLocal * sizeof(int));
    GCSLocal = calloc(totalLocal, sizeof(int));

    if (rank == 0) {
        for (int i = 0; i < totalLocal; i++) {
            tamLocal[i] = tamTotal[i];
            seqLocal[i] = malloc((tamLocal[i] + 1) * sizeof(char));
            strcpy(seqLocal[i], seqTotal[i]);
        }

        offset = totalLocal;
        for (int k = 1; k < size; k++) {
            if (k == size - 1) {
                totaldeSend = sequenciasMPI + resto;
            } else {
                totaldeSend = sequenciasMPI;
            }
            for (int i = 0; i < totaldeSend; i++) {
                idx = offset + i;
                len = tamTotal[idx];

                MPI_Send(&len, 1, MPI_INT, k, 0, MPI_COMM_WORLD);
                MPI_Send(seqTotal[idx], len + 1, MPI_CHAR, k, 1, MPI_COMM_WORLD);
            }
            offset += totaldeSend;
        }
    } else {
        for (int i = 0; i < totalLocal; i++) {
            MPI_Recv(&tamLocal[i], 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            seqLocal[i] = malloc((tamLocal[i] + 1) * sizeof(char));
            MPI_Recv(seqLocal[i], tamLocal[i] + 1, MPI_CHAR, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    #pragma omp parallel for reduction(+:ALocal, CLocal, TLocal, GLocal, GCLocal) schedule(dynamic)
    for (int i = 0; i < totalLocal; i++) {
        int len_seq = tamLocal[i];
        for (int j = 0; j < len_seq; j++) {
            char base = seqLocal[i][j];
            if (base == 'A' || base == 'a') {
                ALocal++;
            } else if (base == 'C' || base == 'c') {
                CLocal++;
                GCLocal++;
                GCSLocal[i]++;
            } else if (base == 'T' || base == 't') {
                TLocal++;
            } else if (base == 'G' || base == 'g') {
                GLocal++;
                GCLocal++;
                GCSLocal[i]++;
            }
        }
    }

    MPI_Reduce(&ALocal,  &ATotal,  1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&CLocal,  &CTotal,  1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&TLocal,  &TTotal,  1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&GLocal,  &GTotal,  1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&GCLocal, &GCTotal, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        for (int i = 0; i < totalLocal; i++) {
            GCSTotal[i] = GCSLocal[i];
        }
        offset = totalLocal;
        for (int l = 1; l < size; l++) {
            if (l == size - 1) {
                totaldeRecv = sequenciasMPI + resto;
            } else {
                totaldeRecv = sequenciasMPI;
            }
            MPI_Recv(&GCSTotal[offset], totaldeRecv, MPI_INT, l, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            offset += totaldeRecv;;
        }
    } else {
        MPI_Send(GCSLocal, totalLocal, MPI_INT, 0, 2, MPI_COMM_WORLD);
    }

    if (rank == 0) {
        printf("\n--- Resumo por entry ---\n");
        for (int i = 0; i < tamanhoTotal; i++) {
            GCLocalAntes = ((double)GCSTotal[i] / tamTotal[i]) * 100.0;

            printf("Seq %d: %s | Tam: %d bases | GC: %d (%.2f%%)\n", i + 1, headersTotal[i], tamTotal[i], GCSTotal[i], GCLocalAntes);
        }

        totalBases = ATotal + CTotal + TTotal + GTotal;
        printf("\n--- Resumo Global ---");
        printf("\nTotal de sequências: %d", tamanhoTotal);
        printf("\nTotal de bases: %d", totalBases);
        printf("\nA: %d\nC: %d\nT: %d\nG: %d", ATotal, CTotal, TTotal, GTotal);
        printf("\nGC Total = %.2f%%\n", ((double)GCTotal / totalBases) * 100.0);

        for (int i = 0; i < tamanhoTotal; i++) {
            free(seqTotal[i]);
        }
        
        free(headersTotal);
        free(seqTotal);
        free(tamTotal);
        free(GCSTotal);
    }

    for (int i = 0; i < totalLocal; i++) {
        free(seqLocal[i]);
    }

    free(seqLocal);
    free(tamLocal);
    free(GCSLocal);

    MPI_Finalize();
    return 0;
}