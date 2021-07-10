/*
 * @(#)FuzzyCMeans.c
 *
 * Copyright (c) 2013 gyaikhom
 *
 */
/*--------------------------------------------------------------
--  #1.
--  Date: Aug, 08, 2019
--  Author: gyaikhom
--  Motivo: Contribute to students
-------------------------------------------------------------
--  #2.
--  Date: Jun, 01, 2021
--  Author: José Mateus e José Lucas
--  Motivo: Trabalho da disciplina de Sistemas Embarcados 2021.1, Professor: Elias Teodoro
-------------------------------------------------------------
**/


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_DATA_POINTS 178
#define MAX_CLUSTER 3
#define MAX_DATA_DIMENSION 13

int num_data_points;
int num_clusters;
int num_dimensions;
double low_high[MAX_DATA_DIMENSION][2];
double degree_of_memb[MAX_DATA_POINTS][MAX_CLUSTER];
double epsilon;
double fuzziness;
double data_point[MAX_DATA_POINTS][MAX_DATA_DIMENSION];
double cluster_centre[MAX_CLUSTER][MAX_DATA_DIMENSION];

int init(char *fname) {
    int i, j, r, rval;
    FILE *f;
    double s;
    if ((f = fopen(fname, "r")) == NULL) {
        printf("Failed to open input file.");
        return -1;
    }
    fscanf(f, "%d %d %d", &num_data_points, &num_clusters, &num_dimensions);
    if (num_clusters > MAX_CLUSTER) {
        printf("Number of clusters should be < %d\n", MAX_CLUSTER);
        goto failure;
    }
    if (num_data_points > MAX_DATA_POINTS) {
        printf("Number of data points should be < %d\n", MAX_DATA_POINTS);
        goto failure;
    }
    if (num_dimensions > MAX_DATA_DIMENSION) {
        printf("Number of dimensions should be >= 1.0 and < %d\n",
                MAX_DATA_DIMENSION);
        goto failure;
    }
    fscanf(f, "%lf", &fuzziness);
    if (fuzziness <= 1.0) {
        printf("Fuzzyness coefficient should be > 1.0\n");
        goto failure;
    }
    fscanf(f, "%lf", &epsilon);
    if (epsilon <= 0.0 || epsilon > 1.0) {
        printf("Termination criterion should be > 0.0 and <= 1.0\n");
        goto failure;
    }

    for (i = 0; i < num_data_points; i++) {
        // printf("\n%i - ", i);
        for (j = 0; j < num_dimensions; j++) {
            fscanf(f, "%lf", &data_point[i][j]);
            // printf("%f ", data_point[i][j]);
            if (data_point[i][j] < low_high[j][0])
                low_high[j][0] = data_point[i][j];
            if (data_point[i][j] > low_high[j][1])
                low_high[j][1] = data_point[i][j];
        }
    }
    printf("\n");
    for (i = 0; i < num_data_points; i++) {
        s = 0.0;
        r = 100;
        for (j = 1; j < num_clusters; j++) {
            rval = rand() % (r + 1);
            r -= rval;
            degree_of_memb[i][j] = rval / 100.0;
            s += degree_of_memb[i][j];
        }
        degree_of_memb[i][0] = 1.0 - s;
    }
    fclose(f);
    return 0;
failure:
    fclose(f);
    exit(1);
}

int calculate_centre_vectors() {
    int i, j, k;
    double numerator, denominator;
    double t[MAX_DATA_POINTS][MAX_CLUSTER];

    for (i = 0; i < num_data_points; i++) {
        for (j = 0; j < num_clusters; j++) {
            t[i][j] = pow(degree_of_memb[i][j], fuzziness);
        }
    }
    for (j = 0; j < num_clusters; j++) {
        for (k = 0; k < num_dimensions; k++) {
            numerator = 0.0;
            denominator = 0.0;
            for (i = 0; i < num_data_points; i++) {
                numerator += t[i][j] * data_point[i][k];
                denominator += t[i][j];
            }
            cluster_centre[j][k] = numerator / denominator;
        }
    }
    return 0;
}

double get_norm(int i, int j) {
    int k;
    double sum = 0.0;
    for (k = 0; k < num_dimensions; k++) {
        sum += pow(data_point[i][k] - cluster_centre[j][k], 2);
    }
    return sqrt(sum);
}

double get_new_value(int i, int j) {
    int k;
    double t, p, sum;
    sum = 0.0;
    p = 2 / (fuzziness - 1);
    for (k = 0; k < num_clusters; k++) {
        t = get_norm(i, j) / get_norm(i, k);
        t = pow(t, p);
        sum += t;
    }
    return 1.0 / sum;
}

double update_degree_of_membership() {
    int i, j;
    double new_uij;
    double max_diff = 0.0, diff;
    for (j = 0; j < num_clusters; j++) {
        for (i = 0; i < num_data_points; i++) {
            new_uij = get_new_value(i, j);
            diff = new_uij - degree_of_memb[i][j];
            if (diff > max_diff)
                max_diff = diff;
            degree_of_memb[i][j] = new_uij;
        }
    }
    return max_diff;
}

int fuzzy_c_means(char *fname) {
    double max_diff;
    init(fname);
    do {
        calculate_centre_vectors();
        max_diff = update_degree_of_membership();
    } while (max_diff > epsilon);
    return 0;
}

void print_data_points(char *fname) {
    int i, j;
    FILE *f;
    if (fname == NULL)
        f = stdout;
    else if ((f = fopen(fname, "w")) == NULL) {
        printf("Cannot create output file.\n");
        exit(1);
    }
    fprintf(f, "Data points:\n");
    for (i = 0; i < num_data_points; i++) {
        printf("Data[%d]: ", i);
        for (j = 0; j < num_dimensions; j++) {
            printf("%.5lf ", data_point[i][j]);
        }
        printf("\n");
    }
    if (fname == NULL)
        fclose(f);
}

void print_membership_matrix(char *fname) {
    int i, j;
    int clusters[3] = {0};
    FILE *f;
    if (fname == NULL)
        f = stdout;
    else if ((f = fopen(fname, "w")) == NULL) {
        printf("Cannot create output file.\n");
        exit(1);
    }
    fprintf(f, "Membership matrix:\n");
    fprintf(f, "Elements\tCluster 01\tCluster 02\tCluster 03\n");
    for (i = 0; i < num_data_points; i++) {
        float bigger = 0;
        char index = 0;
        fprintf(f, "Data[%d]: ", i);
        for (j = 0; j < num_clusters; j++) {
            fprintf(f, "%lf ", degree_of_memb[i][j]);
            if (bigger < degree_of_memb[i][j]) {
                bigger = degree_of_memb[i][j];
                index = j;
            }
        }
        clusters[index]++;
        fprintf(f, "\n");
    }
    printf("------------------------------------------------------------------------\n");
    printf("Grupo 01: %d\nGrupo 02: %d\nGrupo 03: %d\n", clusters[0], clusters[1], clusters[2]);

    if (fname == NULL)
        fclose(f);
}
