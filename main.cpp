#include <conio.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define imax 1000000
#define ERRO 0.00001

//#define ERRO 1e-6
void prodMVet(double **M, double x[], double r[], int TAM)
{
    int i,j;
    for (i=0; i<TAM; i++)
    {
        double soma=0;
        for (j=0; j<TAM; j++)
            soma+=M[i][j]*x[j];
        r[i]=soma;
    }
}

void prodMtransVet(double **M, double x[], double r[], int TAM)
{
    int i,j;
    for (i=0; i<TAM; i++)
    {
        double soma=0;
        for (j=0; j<TAM; j++)
            soma+=M[j][i]*x[j];
        r[i]=soma;
    }
}

void subVet(double a[], double b[], double r[], int TAM)
{
    int i;
    for (i=0; i<TAM; i++) r[i]=a[i]-b[i];
}

void somaVet(double a[], double b[], double r[], int TAM)
{
    int i;
    for (i=0; i<TAM; i++) r[i]=a[i]+b[i];
}

double prodEsc(double a[], double b[], int TAM)
{
    double soma=0;
    int i;
    for (i=0; i<TAM; i++)
        soma+=a[i]*b[i];
    return soma;
}

int gc(double **A, double *b, double *x, int TAM)
{
    double erro=0.00001;
    int i;
    for (i=0; i<TAM; i++) x[i]=0;
    int it=0; // iteracao corrente
    double *aux=(double *)malloc(TAM*sizeof(double));
    prodMVet(A,x,aux,TAM);
    double *r=(double *)malloc(TAM*sizeof(double)); // vetor do res�duo
    subVet(b,aux,r,TAM);
    double *d=(double *)malloc(TAM*sizeof(double));
    for (i=0; i<TAM; i++) d[i]=r[i];
    double sigma_novo=prodEsc(r,r,TAM);
    double sigma0=sigma_novo;
    while (it<imax && sigma_novo>erro*erro*sigma0)
    {
        if(it % 500 == 0)
            printf("it=%d sigma=%lf\n",it,sigma_novo);
        double *q=(double *)malloc(TAM*sizeof(double));
        prodMVet(A,d,q,TAM);
        double Alpha=sigma_novo/prodEsc(d,q,TAM);
        for (i=0; i<TAM; i++) x[i]+=Alpha*d[i];
        if (it%50==0)
        {
            prodMVet(A,x,aux,TAM);
            for (i=0; i<TAM; i++) r[i]=b[i]-aux[i];
        }
        else
            for (i=0; i<TAM; i++) r[i]-=Alpha*q[i];
        double sigma_velho = sigma_novo;
        sigma_novo = prodEsc(r,r,TAM);
        double beta = sigma_novo / sigma_velho;
        for (i=0; i<TAM; i++) d[i]=r[i]+beta*d[i];
        it=it+1;
    }
    return it;
}
double prodEscT(double a[], double b[], int TAM)
{
    double soma=0;
    int i;
    for (i=0; i<TAM; i++)
        soma+=a[TAM-i]*b[i];
    return soma;
}
void mostra(double **A, int TAM)
{
    int i,j;
    for (i=0; i<TAM; i++)
    {
        for (j=0; j<TAM; j++)
            printf("%10lf ",A[i][j]);
        printf("\n");
    }
}

int gauss(double **A, double *b, double *x, int n)
{
    int i,j,k;
    for(j=0; j<n-1; j++)
    {
        for(i=j+1; i<n; i++)
        {
            if (A[j][j]==0) return 0;
            double c=A[i][j]/A[j][j];
            for(k=0; k<n; k++)
                A[i][k]=A[i][k]-c*A[j][k];
            b[i]=b[i]-c*b[j];
        }
    }

    if (A[n-1][n-1]==0) return 0;
    x[n-1]=b[n-1]/A[n-1][n-1];
//    printf("x[%d]=%lf\n",n-1,x[n-1]);
    for(i=n-2; i>=0; i--)
    {
        double sum=0;
        for(j=i+1; j<n; j++)
            sum=sum+A[i][j]*x[j];
        x[i]=(b[i]-sum)/A[i][i];
    }
    return 1;
}
void bicg(double **A, double *b, double *x, int NROW, int NCOL)
{
//Declarações
    int i=0,j=0;

//nlin, ncols = A.shape
    int nLin = NROW;
    int nCols = NCOL;
//print("nlin=",nlin,"ncols=",ncols)
    printf ("Numero Linhas %d | Numero Colunas %d \n",nLin,nCols);
//x = np.zeros((nlin,1))
    for (i=0; i<nCols; i++)
    {
        x[i] = 0;
    }
//r = np.subtract(b,np.matmul(A,x))
    double *aux=(double *)malloc(NCOL*sizeof(double));
    prodMVet(A,x,aux,nCols);
    double *r=(double*)calloc(NCOL,sizeof(double));
    subVet(b,aux,r,nCols);
    //for( i = 0; i<nCols;i++){
    //printf("b[%d]= %lf\n",i,x[i]);
    //r[i] = b[i] - ( *aux);
    //printf("r[%d]= %lf\n",i,r[i]);
    // }

//r2 = r
    double *r2=(double*)calloc(NCOL,sizeof(double));
    for(i=0; i<nCols; i++)
    {
        r2[i] = r[i];
    }
//p = np.zeros((nlin,1))
    double *p=(double*)calloc(NCOL,sizeof(double));
    for(i=0; i<nCols; i++)
    {
        p[i] = 0;
    }
//p2 = np.zeros((nlin,1))
    double *p2=(double*)calloc(NCOL,sizeof(double));
    for(i=0; i<nCols; i++)
    {
        p2[i] = 0;
    }
//rho = 1
    double rho = 1;

//#print(r)
//i = 1
//Declaração
    int it = 1;
    double rho0 = 1;
    double beta= 0;
    double *v=(double *)malloc(NCOL*sizeof(double));
    double alpha = 0.0;
    double erro = 0.0;
    double *aux3=(double *)malloc(NCOL*sizeof(double));
    double aux4 = 0;

//while  i < IMAX:
    while (it < imax)
    {
//rho0 = rho
        rho0 = rho;
//rho = np.dot(np.transpose(r2), r)
        rho = prodEsc(r2,r,nCols);
        //rho = prodEscT(r2,r,nCols);
        //printf("prodESC %lf",prodEsc(r2,r,nCols));
        printf ("RHO %lf - ",rho);
        //for( i = 0; i<nCols;i++){
        //     printf("r[%d]= %lf\n",i,r[i]);
        //}
//beta = rho / rho0
        beta = rho /rho0;
        printf ("BETA %lf -",beta);
//p = np.add(r,np.multiply(beta,p))
        for( i = 0; i<nCols; i++)
        {
            p[i] =r[i]+(beta*p[i]);
            //printf("p[%d]= %lf\n",i,p[i]);
        }
//p2 = np.add(r2,np.multiply(beta,p2))
        for( i = 0; i<nCols; i++)
        {
            p2[i] = r2[i]+ (beta*p2[i]);
            //printf("p2[%d]= %lf\n",i,p2[i]);
        }
//v = np.matmul(A,p)
        prodMVet(A,p,v,nCols);
        //for (i=0;i<nCols;i++){
        //    printf("V[%d]: %f\n" ,i,x[i]);
        // }
//alpha = rho/np.dot(np.transpose(p2), v)
        alpha = rho / prodEsc(p2,v,nCols);
        //alpha = rho / prodEscT(p2,v,nCols);
        //printf("prodESC %lf",prodEsc(p2,v,nCols));
        printf ("ALPHA %lf - ",alpha);

//x = np.add(x,np.multiply(alpha,p))
        for( i = 0; i<nCols; i++)
        {
            aux4 = (alpha*p[i]);
            x[i] = x[i]+ aux4;
            //printf("x[%d]= %lf\n",i,x[i]);
        }
//#print(r)
//erro = np.dot(np.transpose(r), r)
        erro = prodEsc(r,r,nCols);
        //erro = prodEscT(r,r,nCols);
//print(erro)
        printf ("ERRO: %lf\n - ",erro);
//if erro < ERRO * ERRO:
        if (erro < (ERRO * ERRO))
//    break
            break;
//r = np.subtract(r, np.multiply(alpha,v))
        aux4 = 0;
        for(i = 0; i< nCols; i++)
        {
            aux4 = ( alpha * v[i]);
            r[i] = r[i] - aux4;
            //printf("r[%d]=%lf - ",i,r[i]);
        }
//r2 = np.subtract(r2, np.multiply(alpha, np.matmul(np.transpose(A),p2)))
        aux4 = 0;
        prodMtransVet(A,p2,aux3,nCols);
        for(i=0; i<nCols; i++)
        {
            aux4 = ( alpha * (aux3[i]));
            r2[i] = r2[i] - aux4;
            //r2[i] = r2[i] - ( alpha * (*aux3));
            //printf("r2[%d]= %lf\n",i,r2[i];
        }
//i = i + 1
        it+= 1;
        printf("ITERACOES ATUAL: %d\n", it);

    }

//return x, i
    for (i=0; i<nCols; i++)
    {
        printf("Resultado: X[%d] %lf\n",i,x[i]);
    }
    printf("iteracoes: %d",it);
}

void trata_arq(char *arquivo)
{
    FILE *arqin = fopen(arquivo, "rt");
     printf("- %s -\n",arquivo);

    if (!arqin)
    {
        printf("Erro na abertura de %s %d\n",arquivo,strlen(arquivo));
        exit(0);
    }
    char linha[100];
    fgets(linha, 100, arqin);
    printf("%s\n", linha);

    int TOTCRD, PTRCRD, INDCRD, VALCRD, RHSCRD;
    char MXTYPE[3];
    int NROW, NCOL, NNZERO, NELTVL;
    char PTRFMT[20], INDFMT[20], VALFMT[20], RHSFMT[20];

    fgets(linha, 100, arqin);
    sscanf(linha, "%d %d %d %d %d", &TOTCRD, &PTRCRD, &INDCRD, &VALCRD, &RHSCRD);
    printf("TOTCRD = %d - Total number of lines excluding header\n",TOTCRD);
    printf("PTRCRD = %d - Number of lines for pointers\n",PTRCRD);
    printf("INDCRD = %d - Number of lines for variables indices\n",INDCRD);
    printf("VALCRD = %d - Number of lines for numerical values\n",VALCRD);

    fgets(linha, 100, arqin);
    sscanf(linha, "%s %d %d %d %d", MXTYPE, &NROW, &NCOL, &NNZERO, &NELTVL);
    printf("RHSCRD = %d - Number of lines for right-hand sides\n\n",RHSCRD);
    printf("MXTYPE = %s - Matrix type\n",MXTYPE);
    printf("NROW = %d - Number of rows (for variables)\n",NROW);
    printf("NCOL = %d - Number of columns\n",NCOL);
    printf("NNZERO = %d - Number of row indices (elements)\n",NNZERO);
    printf("NELTVL = %d - Number of elemental matrix entries\n\n",NELTVL);

    fgets(linha, 100, arqin);
    printf("%s \n",linha);
    strcpy(PTRFMT,strtok(linha," "));
    printf("%s \n",PTRFMT);
    strcpy(INDFMT,strtok(NULL," "));
    printf("%s\n",INDFMT);
    strcpy(VALFMT,strtok(NULL," "));
    printf("%s\n",VALFMT);
    char *paux;
    if (paux=strtok(NULL," "))
        strcpy(RHSFMT,paux);
    else
        strcpy(RHSFMT,"");
    printf("%s\n",RHSFMT);
    printf("%s \n",PTRFMT);
    printf("%s\n",INDFMT);
    printf("%s\n",VALFMT);
    printf("PTRFMT = %s - Format for pointers\n",PTRFMT);
    printf("INDFMT = %s - Format for indices\n",INDFMT);
    printf("VALFMT = %s - Format for numerical values\n",VALFMT);
    printf("RHSFMT = %s - Format for numerical values of RHS\n\n",RHSFMT);




    int n_linhas1, n_elementos1;
    sscanf(PTRFMT, "(%dI%d)", &n_linhas1, &n_elementos1);


    int n_linhas2, n_elementos2;
    sscanf(INDFMT, "(%dI%d)", &n_linhas2, &n_elementos2);

    printf("Indices: h� %d elementos de %d colunas por linha\n",n_linhas2,n_elementos2);


    int n_elementos_linhas, n_elementos_coluna, tam_elemento;
    if (strstr(VALFMT,"P")){
        sscanf(VALFMT, "(%*dP%dE%d.%d)", &n_elementos_linhas, &n_elementos_coluna, &tam_elemento);

        }
    else
        sscanf(VALFMT, "(%dE%d.%d)", &n_elementos_linhas, &n_elementos_coluna, &tam_elemento);
    printf("Valores: h� %d elementos com tamanho %d e %d decimais\n",
           n_elementos_linhas,n_elementos_coluna,tam_elemento);
    // Inicio dos elementos de cada coluna
    int *inicol = (int *)calloc(NROW+1,sizeof(int));

    // N�mero da linha de cada elemento
    int *numlin = (int *)calloc(NNZERO,sizeof(int));
    double *entries = (double *)calloc(NNZERO,sizeof(double));

    printf("(%dP%dE%d.%d)\n",n_elementos_linhas, n_elementos_coluna, tam_elemento);


    printf("\n\nInicio dos dados de cada coluna\n\n");
    char line[100];
    char dest[6];
    int k = 0;
    char st[100], staux[n_elementos1 + 1];
    char *p;
    int int_aux,i;

    for (i = 0; i < NROW+1; i++)
    {
        if (i % n_linhas1 == 0)
        {
            fgets(line, 100, arqin);
            p = line;
        }
        strncpy(staux, p, n_elementos1);
        staux[n_elementos1] = '\0';
        sscanf(staux, "%d", &int_aux);
        inicol[k] = int_aux-1;
        k++;
        p += n_elementos1;
    }
    k = 0;


    for(int i =0; i<NCOL+1; i++)
    {
        printf("%d ",inicol[i]);
        fflush(stdout);
    }

    for (i = 0; i < NNZERO; i++)
    {

        if (i % n_linhas2 == 0)
        {
            fgets(line, 100, arqin);
            p = line;
        }
        strncpy(staux, p, n_elementos2);

        staux[n_elementos2] = '\0';
        sscanf(staux, "%d", &int_aux);
        numlin[k] = int_aux-1;
        k++;
        p += n_elementos2;
    }
    k = 0;


    /*        printf("\nN�mero da linha de cada elemento\n");
            for(int i = 0;i<NNZERO;i++)
                printf("%d ",numlin[i]);*/

    char staux2[n_elementos_coluna];
    double aux2;
    for (i = 0; i < NNZERO; i++)
    {
        if (i % n_elementos_linhas == 0)
        {
            fgets(line, 100, arqin);
            p = line;
        }

        strncpy(staux2, p, n_elementos_coluna);
        staux[n_elementos_coluna] = '\0';
        sscanf(staux2, "%lf", &aux2);
        entries[i] = aux2;
        p += n_elementos_coluna;
    }
    printf("\nValores\n\n");
    for(int i = 0; i<NNZERO; i++)
        printf("%lf\n",entries[i]);

    printf("\nVou alocar matriz. NROW=%d NCOL=%d\n",NROW,NCOL);
    double **A=(double **)malloc(NROW*sizeof(double *));

    printf("\nVou alocar matriz.\n");
    for(int i = 0; i < NROW; i++)
    {
        A[i] = (double *)calloc(NCOL, sizeof(double));
    }
    printf("\nAloquei matriz\n");
    int a=0;

    int contEntries = 0;
    for (i = 0; i < NCOL; i++)
    {
        //printf("i=%d\n",i);
        int limite = 0;
        if (i != NCOL)
        {
            limite = inicol[i + 1];
        }
        int j,finalRow = i + 1;
        for (j = inicol[i]; j < limite; j++)
        {
            //printf("j=%d\n",i);
            A[numlin[j]][i] = entries[contEntries];
            printf("linha %d coluna %d recebe %lf\n",numlin[j],i,entries[contEntries]);
            //matriz[i][numlin[j]] = entries[contEntries]; //se a matriz � sim�trica, coloca tanto no [i][j] quanto no [j][i]
            contEntries++;
        }
    }

    double *x=(double*)calloc(NCOL,sizeof(double));
    double *b=(double*)calloc(NCOL,sizeof(double));
    for (i=0; i<NCOL; i++)
        b[i]=1.0;

    printf("\n\nVou gerar o arquivo de saida\n\n");
    FILE *saida=fopen("saida.txt","wt");
    for(int i=0; i<NROW; i++)
    {
        for(int j=0; j<NCOL; j++)
        {
            if (A[i][j]!=0)
                fprintf(saida,"X");
            else
                fprintf(saida," ");
            //fprintf(saida,"%16.8f ",A[i][j]);
        }
        fprintf(saida,"\n");
    }
    fclose(saida);

    bicg(A,b,x,NROW,NCOL);
}



main()
{
    char nome[200];
    FILE *arquivos=fopen("arquivo.txt","rt");
    while (fgets(nome,200,arquivos))
    {
        for (int i=0; nome[i]; i++) if (nome[i]=='\n') nome[i]='\0';
        trata_arq(nome);
    }
}
