#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#define PI 3.141592653589793

int main(int argc,char **argv) {
    const double f0 = 38000.0;
    const double c0 = 1480.0;
    const double k = 2.0 * PI * f0 / c0;
    const double d[3] = {1.0, 0.0, 0.0};

    FILE *fp = fopen(argc>1?argv[1]:"sphere-1.905-600.msh", "r");
    if (!fp) {
        perror("Cannot open mesh file");
        return 1;
    }

    char line[1024];
    int n, m;

    // Skip header
    for (int i = 0; i < 4; ++i) fgets(line, sizeof(line), fp);
    fgets(line, sizeof(line), fp);
    sscanf(line, "%d", &n);

    double **v = malloc(n * sizeof(double *));
    for (int i = 0; i < n; ++i) {
        v[i] = malloc(3 * sizeof(double));
        fgets(line, sizeof(line), fp);
        sscanf(line, "%*d %lf %lf %lf", &v[i][0], &v[i][1], &v[i][2]);
    }

    for (int i = 0; i < 2; ++i) fgets(line, sizeof(line), fp);
    fgets(line, sizeof(line), fp);
    sscanf(line, "%d", &m);

    int **e = malloc(m * sizeof(int *));
    for (int i = 0; i < m; ++i) {
        e[i] = malloc(3 * sizeof(int));
        fgets(line, sizeof(line), fp);
        sscanf(line, "%*d %*d %*d %*d %*d %d %d %d",
               &e[i][0], &e[i][1], &e[i][2]);
        e[i][0]--; e[i][1]--; e[i][2]--; // Convert to 0-based indexing
    }

    fclose(fp);

    // Triangle centers
    double **x = malloc(m * sizeof(double *));
    #pragma omp parallel for
    for (int i = 0; i < m; ++i) {
        x[i] = calloc(3, sizeof(double));
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                x[i][k] += v[e[i][j]][k] / 3.0;
    }
    
    // Helmholtz matrix
    double complex **S = malloc(m * sizeof(double complex *));
    #pragma omp parallel for
    for (int i = 0; i < m; ++i) {
        S[i] = malloc(m * sizeof(double complex));
        for (int j = 0; j < m; ++j) {
            double r2 = 0.0;
            for (int d = 0; d < 3; ++d)
                r2 += pow(x[i][d] - x[j][d], 2);
            double r = sqrt(r2);
            S[i][j] = (r != 0) ? cexp(I * k * r) / r / (4.0 * PI)
                               : I * k / (4.0 * PI);
        }
    }

    // Incident field
    double complex *f = malloc(m * sizeof(double complex));
    #pragma omp parallel for
    for (int i = 0; i < m; ++i) {
        /*
        double dot = 0.0;
        for (int j = 0; j < 3; ++j)
            dot += x[i][j] * d[j];
        f[i] = cexp(I * k * dot);
        */
        f[i] = cexp(I * k * (x[i][0]*d[0]+x[i][1]*d[1]+x[i][2]*d[2]) );
    }

    // Negate f
    double complex *f_neg = malloc(m * sizeof(double complex));
    for (int i = 0; i < m; ++i)
        f_neg[i] = -f[i];

    // Gaussian elimination
    double complex *phi = calloc(m, sizeof(double complex));
    for (int k = 0; k < m; ++k) {
        int max = k;
        for (int i = k + 1; i < m; ++i)
            if (cabs(S[i][k]) > cabs(S[max][k])) max = i;
        if (max != k) {
            double complex *tmp = S[k]; S[k] = S[max]; S[max] = tmp;
            double complex tmp_f = f_neg[k]; f_neg[k] = f_neg[max]; f_neg[max] = tmp_f;
        }
        for (int i = k + 1; i < m; ++i) {
            double complex factor = S[i][k] / S[k][k];
            for (int j = k; j < m; ++j)
                S[i][j] -= factor * S[k][j];
            f_neg[i] -= factor * f_neg[k];
        }
    }
    for (int i = m - 1; i >= 0; --i) {
        phi[i] = f_neg[i];
        for (int j = i + 1; j < m; ++j)
            phi[i] -= S[i][j] * phi[j];
        phi[i] /= S[i][i];
    }

    // Far field
    #pragma omp parallel for
    for (int i = 0; i < 360; ++i) {
        double theta = PI * i / 180.0;
        double r[3] = {cos(theta), sin(theta), 0.0};
        double complex sum = 0.0 + 0.0 * I;
        for (int j = 0; j < m; ++j) {
            double dot = 0.0;
            for (int d_idx = 0; d_idx < 3; ++d_idx)
                dot += x[j][d_idx] * r[d_idx];
            sum += cexp(-I * k * dot) / (4.0 * PI) * phi[j];
        }
        printf("%d\t%g\n", i, cabs(sum));
    }

    // Cleanup
    for (int i = 0; i < n; ++i) free(v[i]);
    for (int i = 0; i < m; ++i) {
        free(e[i]); free(x[i]); free(S[i]);
    }
    free(v); free(e); free(x); free(S);
    free(f); free(f_neg); free(phi);

    return 0;
}
