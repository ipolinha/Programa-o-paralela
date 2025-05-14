#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#define NX 32
#define NY 32
#define NZ 32
#define NSTEPS 10000
#define DX 1.0
#define DY 1.0
#define DZ 1.0
#define DT 0.005
#define NU 0.05
#define SAVE_EVERY 1000

int main() {
    double u[NX][NY][NZ], v[NX][NY][NZ], w[NX][NY][NZ];
    double u_new[NX][NY][NZ], v_new[NX][NY][NZ], w_new[NX][NY][NZ];

    // Inicialização
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                u[i][j][k] = 0.0;
                v[i][j][k] = 0.0;
                w[i][j][k] = 0.0;
            }
        }
    }

    // Perturbação inicial no centro
    //int cx = NX / 2, cy = NY / 2, cz = NZ / 2;
    //u[cx][cy][cz] = 1.0;

    struct timeval start, end;
    gettimeofday(&start, NULL);  // Início da medição

    FILE* fp = fopen("velocidade_u_sempertubacao.csv", "w");
    if (!fp) {
        perror("Erro ao abrir arquivo de saída");
        exit(1);
    }

    // Cabeçalho CSV
    fprintf(fp, "tempo,i,j,k,u\n");

    //Evolução no tempo
    for (int n = 0; n <= NSTEPS; n++) {
        // Atualização com difusão
        for (int i = 1; i < NX - 1; i++) {
            for (int j = 1; j < NY - 1; j++) {
                for (int k = 1; k < NZ - 1; k++) {
                    //Calcula o Laplaciano de cada componente com diferenças finitas centradas e aplica a difusão u(t+1) = u(t) + ν * Δt * ∇²u
                    double lap_u = (u[i+1][j][k] - 2*u[i][j][k] + u[i-1][j][k]) / (DX*DX)
                                 + (u[i][j+1][k] - 2*u[i][j][k] + u[i][j-1][k]) / (DY*DY)
                                 + (u[i][j][k+1] - 2*u[i][j][k] + u[i][j][k-1]) / (DZ*DZ);

                    double lap_v = (v[i+1][j][k] - 2*v[i][j][k] + v[i-1][j][k]) / (DX*DX)
                                 + (v[i][j+1][k] - 2*v[i][j][k] + v[i][j-1][k]) / (DY*DY)
                                 + (v[i][j][k+1] - 2*v[i][j][k] + v[i][j][k-1]) / (DZ*DZ);

                    double lap_w = (w[i+1][j][k] - 2*w[i][j][k] + w[i-1][j][k]) / (DX*DX)
                                 + (w[i][j+1][k] - 2*w[i][j][k] + w[i][j-1][k]) / (DY*DY)
                                 + (w[i][j][k+1] - 2*w[i][j][k] + w[i][j][k-1]) / (DZ*DZ);

                    u_new[i][j][k] = u[i][j][k] + NU * DT * lap_u;
                    v_new[i][j][k] = v[i][j][k] + NU * DT * lap_v;
                    w_new[i][j][k] = w[i][j][k] + NU * DT * lap_w;
                }
            }
        }

        // Atualiza os valores da próxima etapa para a etapa atual.
        for (int i = 1; i < NX - 1; i++) {
            for (int j = 1; j < NY - 1; j++) {
                for (int k = 1; k < NZ - 1; k++) {
                    u[i][j][k] = u_new[i][j][k];
                    v[i][j][k] = v_new[i][j][k];
                    w[i][j][k] = w_new[i][j][k];
                }
            }
        }

        // Salvar todo o volume
        //A cada SAVE_EVERY passos, salva todos os valores de u para o tempo atual no arquivo.
        //Grava cada linha com: tempo, i, j, k, valor_u.
        if (n % SAVE_EVERY == 0) {
            for (int i = 0; i < NX; i++) {
                for (int j = 0; j < NY; j++) {
                    for (int k = 0; k < NZ; k++) {
                        fprintf(fp, "%d,%d,%d,%d,%.5f\n", n, i, j, k, u[i][j][k]);
                    }
                }
            }
            printf("Salvo passo %d\n", n);
        }
    }

    fclose(fp);
    gettimeofday(&end, NULL);  // Fim da medição

    double tempo = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1e6;
    printf("Tempo de execução: %.6f segundos\n", tempo);

    return 0;
}
