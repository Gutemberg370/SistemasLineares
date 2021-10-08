#include <stdio.h>
#include <stdlib.h>
#include "Métodos.h"

int main()
{
    int ordem,i,j,k, Tipo;
    double **Matriz;
    double *resposta;
    char Tipos_De_Solucao[3][25] = {"incompativel","compativel determinado","compativel indeterminado"};

    printf("Digite a ordem da Matriz de coeficientes:\n");
    scanf("%d", &ordem);
    while(ordem <= 0){
        printf("Digite uma ordem diferente de zero:\n");
        scanf("%d", &ordem);
    }

    resposta = malloc(sizeof(double)*ordem);
    Matriz = alocaMatriz(ordem, ordem+1);

    if(Matriz == NULL){
        printf("Memoria insuficiente.");
        return 1;
    }

    for(i = 0; i < ordem; i++){
            resposta[i] = 0;
        for(j = 0; j < ordem+1; j++){
                printf("Digite o termo M[%d][%d] da Matriz Aumentada:\n",i+1,j+1);
                scanf("%lf", &Matriz[i][j]);
        }
    }

    printf("Matriz Aumentada:\n\n");
    for(i = 0; i < ordem; i++){
        for(j = 0; j < ordem+1; j++){
                printf("%10.3lf ", Matriz[i][j]);
        }
        printf("\n\n");
    }

    Tipo = Metodo_De_Gauss(ordem,Matriz,resposta);
    printf("A solu%c%co %c do tipo %s.\n\n",135,198,130,Tipos_De_Solucao[Tipo]);

    if(Tipo){
        for(k =0; k<ordem; k++)
        {
            printf("x%d = %.4lf; ",k+1, resposta[k]);
        }
        printf("\n\n");
    }

    return 0;
}


