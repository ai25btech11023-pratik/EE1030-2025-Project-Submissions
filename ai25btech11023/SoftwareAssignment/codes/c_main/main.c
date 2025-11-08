#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../c_libs/stb_image.h"
#include "../c_libs/stb_image_write.h"


void mat_mul(int m, int n, const double *A, double *C) {

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            C[i*n + j] = 0.0;
        }
    }

    for (int i = 0; i < n; i++) {          
        for (int j = 0; j < n; j++) {      
            double sum = 0.0;
            for (int k = 0; k < m; k++) {  
                sum += A[k*n + i] * A[k*n + j];
            }
            C[i*n + j] = sum;
        }
    }
}

double normvec(double *a, int n)
{
    double norm = 0.0;
    for (int i = 0; i < n; i++)
        norm += a[i] * a[i];
    return sqrt(norm);
}

void makeUnit(double *a, int n)
{
    double norm = normvec(a, n);
    if (norm == 0)
        return;
    for (int i = 0; i < n; i++)
        a[i] /= norm;
}


//b <--- A @ b
void matvec(int m, int n, double *A, double *x, double *y)
{
    for (int i = 0; i < m; i++)
        y[i] = 0;
    for (int i = 0; i < m; i++)
    {
        double sum = 0.0;
        double *row = A + i * n;
        for (int j = 0; j < n; j++)
            sum += row[j] * x[j];
        y[i] = sum;
    }
}


// b_k = A @ b_k-1
void powerIteration(int n, double *A, double *b, double *sigma, int iterate, double *b_o, double *Ab){
    for(int i=0;i<n;i++) b_o[i]=1;      //initialising b_o as {1,1,1....}
    while(iterate--){
        matvec(n,n,A,b_o,b);
        makeUnit(b,n);
        for(int i=0;i<n;i++) b_o[i]=b[i];
    }
    matvec(n,n,A,b,Ab);
    *sigma = normvec(Ab,n);
}

//A<--- A - sigma *(b @ b^\top)

void subContrib(int n, double *A, double *b, double sigma){
    for(int i=0;i<n;i++){
        double *row = A + i*n;
        for(int j=0;j<n;j++){
            row[j]-= sigma * b[j]*b[i];
        }
    }
}

int main(void){
    char *infilename = (char*)malloc(100*sizeof(char));
    char *outfilename = (char*)malloc(100*sizeof(char));
    printf("Enter input filename: ");
    char infile[100] = "../../figs/";
    char outfile[100] = "../../figs/";
    scanf("%s", infilename);
    strcat(infile, infilename);
    printf("Enter output filename: ");
    scanf("%s", outfilename);
    strcat(outfile, outfilename);
    int m, n, pixint;
    unsigned char *img = stbi_load(infile, &n, &m, &pixint, 1);
    if (!img)
    {
        printf("Failed to load %s\n", infile);
        return 1;
    }
    int k;
    printf("Enter k: ");
    scanf("%d", &k);
    int iterate = 200;

    double *A = (double *)malloc(m * n * sizeof(double));
    double *A_final = (double *)calloc(m * n, sizeof(double));
    double *AtA = (double *)calloc(n * n, sizeof(double));
    double *u = (double *)malloc(m * sizeof(double));
    double *b = (double *)malloc(n * sizeof(double));

    for(int i=0;i<m*n;i++){
        A[i]=img[i];
    }

    mat_mul(m,n,A,AtA);
    double *b_o = (double*)malloc(n* sizeof(double));
    double *Ab = (double *)malloc(n*sizeof(double));
    while(k--){
        double sigma;
        powerIteration(n,AtA,b,&sigma,iterate,b_o,Ab);
	double sqrt_sigma = sqrt(sigma);
        matvec(m,n,A,b,u);
        for (int i = 0; i < m; i++)
        if(sigma!=0) u[i] = u[i] / sqrt_sigma;

        //filling A_final matrix
        for (int i = 0; i < m; i++)
        {
            double ui = u[i];
            double *rowC = A_final + i * n;
            for (int j = 0; j < n; j++)
            {
                rowC[j] += sqrt_sigma * ui * b[j];
            }
        }

        subContrib(n,AtA,b,sigma);
    }
    unsigned char *final = malloc(m * n);
for (int i = 0; i < m * n; i++) {
    double v = A_final[i];
    if (v < 0) v = 0;
    if (v > 255) v = 255;
    final[i] = (unsigned char)round(v);
}
    double errorsum =0;
    for(int i=0;i<m*n;i++){
        errorsum+= pow(A[i]-A_final[i],2);
    }
    errorsum = sqrt(errorsum);
    printf("%lf\n", errorsum);
    stbi_write_jpg(outfile, n, m, 1, final, 50);
    free(A);
    free(A_final);
    free(AtA);
    free(u);
    free(b);
    free(infilename);
    free(outfilename);
    free(final);
    free(b_o);
    free(Ab);
    stbi_image_free(img);
    return 0;
}
