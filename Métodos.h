#include <math.h>
#define EPSILON 0.000000001

////////////////////////////////////////// FUNÇÕES AUXILIARES ////////////////////////////////////////


// Alocação de memória manual para a matriz aumentada
double **alocaMatriz(int l, int c)
{
    // Se houver memória suficiente, aloca uma matriz de double com l linhas e c colunas e devolve
    // um ponteiro para essa matriz; caso contrário devolve um ponteiro nulo.
    double **m;
    int i, j;

    m = malloc(sizeof(double *) * l);
    if (m == NULL)
    {
        return NULL;
    } //memoria insuficiente
    for (i = 0; i < l; i++)
    {
        m[i] = malloc(sizeof(double) * c);
        if (m[i] == NULL)
        {
            for (j = 0; j < i; j++)
            {
                free(m[j]);
            }
            free(m);
            return NULL;
        }
    }
    return m;
}




// Implementação do método de substituição retroativa
int substituicao_retroativa(int ordem, double **Matriz_Aumentada,double *solucao){
    // Recebe como parâmetros a ordem da matriz de coeficientes do sistema linear, sua respectiva Matriz Aumentada na forma
    // triangular superior e um vetor de comprimento igual ao número de variáveis a ser calculado, para compor os resultados obitidos.
    // A função retorna um inteiro identificando o tipo de sistema linear: 0 - não compatível, 1 - compatível determinado
    // e 2 - compatível indeterminado e preenche o vetor solução com a solução do sistema, caso ela exista. No caso de um
    // sistema compatível indeterminado, ele preenche com uma das soluções possíveis.
    int i,j, tipo_de_solucao = 1;
    double Termo_independente;

    for(i = ordem-1; i >= 0 ; i--)
    {
        Termo_independente = Matriz_Aumentada[i][ordem];
        for(j = i+1; j < ordem; j++)
        {
            Termo_independente -= (Matriz_Aumentada[i][j]*solucao[j]);
        }
        if(Matriz_Aumentada[i][i] == 0){
                if(fabs(Termo_independente) < EPSILON){
                  solucao[i] = 0;
                  tipo_de_solucao = 2;
                }
                else{
                    tipo_de_solucao = 0;
                    return tipo_de_solucao;
                }

        }
        else{
            solucao[i] = Termo_independente/ Matriz_Aumentada[i][i];
        }
    }
    return tipo_de_solucao;

}

// Implementação do método de substituição progressiva
int substituicao_progressiva(int ordem, double **Matriz_Aumentada, double *solucao){
    // Recebe como parâmetros a ordem da matriz de coeficientes do sistema linear, sua respectiva Matriz Aumentada na forma
    // triangular inferior e um vetor de comprimento igual ao número de variáveis a ser calculado, para compor os resultados obitidos.
    // A função retorna um inteiro identificando o tipo de sistema linear: 0 - não compatível, 1 - compatível determinado
    // e 2 - compatível indeterminado e preenche o vetor solução com a solução do sistema, caso ela exista. No caso de um
    // sistema compatível indeterminado, ele preenche com uma das soluções possíveis.
    int i,j, tipo_de_solucao = 1;
    double Termo_independente;

    for(i = 0; i < ordem ; i++)
    {
        Termo_independente = Matriz_Aumentada[i][ordem];
        for(j = 0; j < i; j++)
        {
            Termo_independente -= (Matriz_Aumentada[i][j]*solucao[j]);
        }
        if(Matriz_Aumentada[i][i] == 0){
                if(fabs(Termo_independente) < EPSILON){
                  solucao[i] = 0;
                  tipo_de_solucao = 2;
                }
                else{
                    tipo_de_solucao = 0;
                    return tipo_de_solucao;
                }

        }
        else{
            solucao[i] = Termo_independente/ Matriz_Aumentada[i][i];
        }
    }
    return tipo_de_solucao;

}

// Implementação do método de substituição diagonal
int substituicao_diagonal(int ordem, double **Matriz_Aumentada, double *solucao){
    // Recebe como parâmetros a ordem da matriz de coeficientes do sistema linear, sua respectiva Matriz Aumentada na forma
    // triangular inferior e um vetor de comprimento igual ao número de variáveis a ser calculado, para compor os resultados obitidos.
    // A função retorna um inteiro identificando o tipo de sistema linear: 0 - não compatível, 1 - compatível determinado
    // e 2 - compatível indeterminado e preenche o vetor solução com a solução do sistema, caso ela exista. No caso de um
    // sistema compatível indeterminado, ele preenche com uma das soluções possíveis.
    int i, tipo_de_solucao = 1;
    double Termo_independente;

    for(i = 0; i < ordem ; i++)
    {
        Termo_independente = Matriz_Aumentada[i][ordem];
        if(Matriz_Aumentada[i][i] == 0){
                if(fabs(Termo_independente) < EPSILON){
                  solucao[i] = 0;
                  tipo_de_solucao = 2;
                }
                else{
                    tipo_de_solucao = 0;
                    return tipo_de_solucao;
                }

        }
        else{
            solucao[i] = Termo_independente/ Matriz_Aumentada[i][i];
        }
    }
    return tipo_de_solucao;

}

////////////////////////////////////// Métodos Diretos ////////////////////////////////////

// Implementação do método de Gauss
int Metodo_De_Gauss(int ordem, double **Matriz_Aumentada,double *solucao){
    int i,j,k;
    double Mult, *auxiliar;
    double **Matriz_Aumentada_Controle;

    Matriz_Aumentada_Controle = alocaMatriz(ordem, ordem+1);

    for(i = 0; i < ordem; i++){
        for(j = 0; j < ordem+1; j++){
            Matriz_Aumentada_Controle[i][j] = Matriz_Aumentada[i][j];
        }
    }

    for(i = 0; i < ordem-1; i++)
    {
        j = i;
        if(Matriz_Aumentada_Controle[j][i] == 0){
            while(j < ordem && Matriz_Aumentada_Controle[j][i] == 0){
                j++;
            }
            if(j < ordem){
                auxiliar = Matriz_Aumentada_Controle[i];
                Matriz_Aumentada_Controle[i] = Matriz_Aumentada_Controle[j];
                Matriz_Aumentada_Controle[j] = auxiliar;
            }
        }
        if(Matriz_Aumentada_Controle[i][i] != 0){
            for(j = i+1; j < ordem; j++)
            {
                Mult = -(Matriz_Aumentada_Controle[j][i]/Matriz_Aumentada_Controle[i][i]);
                Matriz_Aumentada_Controle[j][i] = 0;

                for(k = i+1; k < ordem+1; k++){
                    Matriz_Aumentada_Controle[j][k] = Mult* Matriz_Aumentada_Controle[i][k] + Matriz_Aumentada_Controle[j][k];
                }
            }
        }
        }

        printf("Matriz Escalonada:\n\n");
        for(i = 0; i < ordem; i++){
            for(j = 0; j < ordem+1; j++){
                printf("%10.3lf ", Matriz_Aumentada_Controle[i][j]);
        }
        printf("\n\n");
        }
        return substituicao_retroativa(ordem,Matriz_Aumentada_Controle,solucao);
}

// Implementação do método de Jordan
int Metodo_De_Jordan(int ordem, double **Matriz_Aumentada,double *solucao){
    int i,j,k, t = 0, resposta;
    int **Matriz_Trocas;
    double Mult, auxiliar;
    double **Matriz_Aumentada_Controle;


    Matriz_Trocas = alocaMatriz(ordem,2);

    Matriz_Aumentada_Controle = alocaMatriz(ordem, ordem+1);

    for(i = 0; i < ordem; i++){
        for(j = 0; j < ordem+1; j++){
            Matriz_Aumentada_Controle[i][j] = Matriz_Aumentada[i][j];
        }
    }

    for(i = 0; i < ordem; i++)
    {
        j = i;
        if(Matriz_Aumentada_Controle[i][j] == 0){
            while(j < ordem && Matriz_Aumentada_Controle[i][j] == 0){
                j++;
            }
            if(j < ordem){
                for(k = 0; k < ordem; k++){
                auxiliar = Matriz_Aumentada_Controle[k][i];
                Matriz_Aumentada_Controle[k][i] = Matriz_Aumentada_Controle[k][j];
                Matriz_Aumentada_Controle[k][j] = auxiliar;
                Matriz_Trocas[t][0] = i;
                Matriz_Trocas[t][1] = j;
                t++;
                }
            }
            else{
                for(k = 0; k < ordem; k++){
                Matriz_Aumentada_Controle[k][i] = 0;
                }
            }
        }
        if(Matriz_Aumentada_Controle[i][i] != 0){
            for(j = 0; j < ordem; j++)
            {
                if(j != i){
                Mult = -(Matriz_Aumentada_Controle[j][i]/Matriz_Aumentada_Controle[i][i]);
                Matriz_Aumentada_Controle[j][i] = 0;

                for(k = i+1; k < ordem+1; k++){
                    Matriz_Aumentada_Controle[j][k] = Mult* Matriz_Aumentada_Controle[i][k] + Matriz_Aumentada_Controle[j][k];
                }
                }
            }
        }

        }

        printf("Matriz Diagonal:\n\n");
        for(i = 0; i < ordem; i++){
            for(j = 0; j < ordem+1; j++){
                printf("%10.3lf ", Matriz_Aumentada_Controle[i][j]);
        }
        printf("\n\n");
        }
        resposta = substituicao_diagonal(ordem,Matriz_Aumentada_Controle,solucao);
        i = 0;
        while(i < t){
            auxiliar = solucao[Matriz_Trocas[i][0]];
            solucao[Matriz_Trocas[i][0]] = solucao[Matriz_Trocas[i][1]];
            solucao[Matriz_Trocas[i][1]] = auxiliar;
            i++;
        }
        return resposta;
}

// Implementação do método da Pivotação Completa
int Metodo_Da_Pivotacao_Completa(int ordem, double **Matriz_Aumentada,double *solucao){
    int i,j,k,l, t = 0,maior,resultado, linha_maior, coluna_maior;
    int Controle_linhas[ordem], Controle_colunas[ordem], Controle_Matriz[ordem];
    double Mult;
    double **Matriz_Aumentada_Controle, **Matriz_Aumentada_Arrumada, *solucao_temporaria;

    solucao_temporaria = malloc(sizeof(double)*ordem);

    Matriz_Aumentada_Controle = alocaMatriz(ordem, ordem+1);

    Matriz_Aumentada_Arrumada = alocaMatriz(ordem,ordem+1);

    for(i = 0; i < ordem; i++){
        for(j = 0; j < ordem+1; j++){
            Matriz_Aumentada_Controle[i][j] = Matriz_Aumentada[i][j];
        }
    }



    for(i = 0; i < ordem-1; i++)
    {
        maior =  fabs(Matriz_Aumentada_Controle[0][0]);
        linha_maior = 0;
        coluna_maior = 0;
        for(k = 0; k < ordem; k++){
                for(l = 0; l < ordem; l++){
                    if(fabs(Matriz_Aumentada_Controle[k][l]) > maior && Controle_Matriz[k] != 1){
                        maior = fabs(Matriz_Aumentada_Controle[k][l]);
                        linha_maior = k;
                        coluna_maior = l;
                    }
            }
        }

        if(maior != 0){
            for(j = 0; j < ordem; j++)
            {
                if(j != linha_maior && Controle_Matriz[j] != 1){
                Mult = -(Matriz_Aumentada_Controle[j][coluna_maior]/Matriz_Aumentada_Controle[linha_maior][coluna_maior]);
                Matriz_Aumentada_Controle[j][coluna_maior] = 0;

                for(k = 0; k < ordem+1; k++){
                        if(k != coluna_maior){
                            Matriz_Aumentada_Controle[j][k] = Mult* Matriz_Aumentada_Controle[linha_maior][k] + Matriz_Aumentada_Controle[j][k];
                        }
                }
                }
            }
            Controle_linhas[t] = linha_maior;
            Controle_colunas[t] = coluna_maior;
            Controle_Matriz[linha_maior] = 1;
            t++;
        }

        }

        for(i = 0; i < ordem; i++){
            if(Controle_Matriz[i] != 1){
                Controle_linhas[t] = i;
                break;
            }
        }
        j = 0;
        k = 0;
        for(i = 0; i < ordem-1; i++){
            j += Controle_colunas[i];
            k += i;
        }
        k += ordem - 1;
        Controle_colunas[t] = k - j;

        printf("%d %d ", j, k);

        printf("Matriz Triangular Superior Desorganizada:\n\n");
        for(i = 0; i < ordem; i++){
            for(j = 0; j < ordem+1; j++){
                printf("%10.3lf ", Matriz_Aumentada_Controle[i][j]);
        }
        printf("\n\n");
        }

        for(i = 0; i < ordem; i++){
            for(j = 0; j < ordem+1; j++){
                if(j != ordem){
                 Matriz_Aumentada_Arrumada[i][j] = Matriz_Aumentada_Controle[Controle_linhas[i]][Controle_colunas[j]];
                }
                else{
                 Matriz_Aumentada_Arrumada[i][j] = Matriz_Aumentada_Controle[Controle_linhas[i]][j];
                }
            }
        }

        printf("Matriz Triangular Superior Organizada:\n\n");
        for(i = 0; i < ordem; i++){
            for(j = 0; j < ordem+1; j++){
                printf("%10.3lf ", Matriz_Aumentada_Arrumada[i][j]);
        }
        printf("\n\n");
        }

        resultado = substituicao_retroativa(ordem,Matriz_Aumentada_Arrumada,solucao_temporaria);
        for(i = 0; i < ordem; i++){
            solucao[Controle_colunas[i]] = solucao_temporaria[i];
        }

        return resultado;

}


////////////////////////////////////// Métodos Iterativos ////////////////////////////////////


// Implementação do método de Jacob
int Metodo_de_Jacob(int ordem, double **Matriz_Aumentada,double *solucao, double Margem_erro, int iteracoes){

    int i,j, convergente_linha = 1, convergente_coluna = 1, tipo_de_solucao = 0, iteracoes_dadas = 0;
    double *Nova_solucao;
    double Termo_independente,soma_linha, soma_coluna,controle, maior_residuo;

    Nova_solucao = malloc(sizeof(double)*ordem);

    for(i = 0; i < ordem; i++){
            controle = fabs(Matriz_Aumentada[i][i]);
        for(j = 0; j < ordem; j++){
            if(j != i){
            soma_coluna += fabs(Matriz_Aumentada[j][i]);
            soma_linha += fabs(Matriz_Aumentada[i][j]);
            }
        }
        if(controle < soma_linha){
            //printf("%lf %lf %lf",controle, soma_coluna, soma_linha);
            convergente_linha--;
        }
        if(controle < soma_coluna){
            //printf("%lf %lf %lf",controle, soma_coluna, soma_linha);
            convergente_coluna--;
        }
        soma_linha = 0;
        soma_coluna = 0;
    }

    if(convergente_linha == 1 || convergente_coluna == 1){
        for(i = 0; i < ordem ; i++)
        {
            Termo_independente = Matriz_Aumentada[i][ordem];
            for(j = 0; j < ordem; j++)
            {
                if(j != i){
                Termo_independente -= (Matriz_Aumentada[i][j]*solucao[j]);
                }
            }
            if(Matriz_Aumentada[i][i] == 0){
                Nova_solucao[i] = 0;
                tipo_de_solucao = 2;
            }
            else{
                Nova_solucao[i] = Termo_independente/ Matriz_Aumentada[i][i];
            }
        }

        maior_residuo = fabs(Nova_solucao[0] - solucao[0]);
        for(i = 1; i < ordem ; i++){
            if(fabs(Nova_solucao[i] - solucao[i]) > maior_residuo){
                maior_residuo = fabs(Nova_solucao[i] - solucao[i]);
            }
        }
        iteracoes_dadas++;

        while(maior_residuo > Margem_erro && iteracoes_dadas < iteracoes){

            for(i = 0; i < ordem ; i++)
            {
                solucao[i] = Nova_solucao[i];
                printf("x%d = %lf /",i, solucao[i]);
            }
            printf("\n\n");

            for(i = 0; i < ordem ; i++)
            {
                Termo_independente = Matriz_Aumentada[i][ordem];
                for(j = 0; j < ordem; j++)
                {
                    if(j != i){
                    Termo_independente -= (Matriz_Aumentada[i][j]*solucao[j]);
                    }
                }
                if(Matriz_Aumentada[i][i] == 0){
                    Nova_solucao[i] = 0;
                    tipo_de_solucao = 2;
                }
                else{
                    Nova_solucao[i] = Termo_independente/ Matriz_Aumentada[i][i];
                }
            }

            maior_residuo = fabs(Nova_solucao[0] - solucao[0]);
            for(i = 1; i < ordem ; i++){
                if(fabs(Nova_solucao[i] - solucao[i]) > maior_residuo){
                    maior_residuo = fabs(Nova_solucao[i] - solucao[i]);
                }
            }
            iteracoes_dadas++;
        }
        tipo_de_solucao = 1;

    }
    return tipo_de_solucao;

}


// Implementação do método de Gauss-Seidel
int Metodo_De_Gauss_Seidel(int ordem, double **Matriz_Aumentada,double *solucao, double Margem_erro, int iteracoes){

    int i,j, convergente_linha = 1, convergente_coluna = 1, tipo_de_solucao = 0, iteracoes_dadas = 0;
    double *Nova_solucao;
    double Termo_independente,soma_linha, soma_coluna,controle, maior_residuo;

    Nova_solucao = malloc(sizeof(double)*ordem);

    for(i = 0; i < ordem; i++){
            controle = fabs(Matriz_Aumentada[i][i]);
            Nova_solucao[i] = solucao[i];
        for(j = 0; j < ordem; j++){
            if(j != i){
            soma_coluna += fabs(Matriz_Aumentada[j][i]);
            soma_linha += fabs(Matriz_Aumentada[i][j]);
            }
        }
        if(controle < soma_linha){
            //printf("%lf %lf %lf",controle, soma_coluna, soma_linha);
            convergente_linha--;
        }
        if(controle < soma_coluna){
            //printf("%lf %lf %lf",controle, soma_coluna, soma_linha);
            convergente_coluna--;
        }
        soma_linha = 0;
        soma_coluna = 0;
    }

    if(convergente_linha == 1 || convergente_coluna == 1){
        for(i = 0; i < ordem ; i++)
        {
            Termo_independente = Matriz_Aumentada[i][ordem];
            for(j = 0; j < ordem; j++)
            {
                if(j != i){
                Termo_independente -= (Matriz_Aumentada[i][j]*Nova_solucao[j]);
                }
            }
            if(Matriz_Aumentada[i][i] == 0){
                Nova_solucao[i] = 0;
                tipo_de_solucao = 2;
            }
            else{
                Nova_solucao[i] = Termo_independente/ Matriz_Aumentada[i][i];
            }
        }

        maior_residuo = fabs(Nova_solucao[0] - solucao[0]);
        for(i = 1; i < ordem ; i++){
            if(fabs(Nova_solucao[i] - solucao[i]) > maior_residuo){
                maior_residuo = fabs(Nova_solucao[i] - solucao[i]);
            }
        }
        iteracoes_dadas++;

        while(maior_residuo > Margem_erro && iteracoes_dadas < iteracoes){

            for(i = 0; i < ordem ; i++)
            {
                solucao[i] = Nova_solucao[i];
                printf("x%d = %lf /",i, solucao[i]);
            }
            printf("\n\n");

            for(i = 0; i < ordem ; i++)
            {
                Termo_independente = Matriz_Aumentada[i][ordem];
                for(j = 0; j < ordem; j++)
                {
                    if(j != i){
                    Termo_independente -= (Matriz_Aumentada[i][j]*Nova_solucao[j]);
                    }
                }
                if(Matriz_Aumentada[i][i] == 0){
                    Nova_solucao[i] = 0;
                    tipo_de_solucao = 2;
                }
                else{
                    Nova_solucao[i] = Termo_independente/ Matriz_Aumentada[i][i];
                }
            }

            maior_residuo = fabs(Nova_solucao[0] - solucao[0]);
            for(i = 1; i < ordem ; i++){
                if(fabs(Nova_solucao[i] - solucao[i]) > maior_residuo){
                    maior_residuo = fabs(Nova_solucao[i] - solucao[i]);
                }
            }
            iteracoes_dadas++;
        }
        tipo_de_solucao = 1;

    }
    return tipo_de_solucao;

}




