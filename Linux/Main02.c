#include "FuzzyCMeans02.c"

int main(int argc, char **argv) {
    printf("------------------------------------------------------------------------\n");
    fuzzy_c_means();
    printf("Number of data points: %d\n", num_data_points);
    printf("Number of clusters: %d\n", num_clusters);
    printf("Number of data-point dimensions: %d\n", num_dimensions);
    printf("Accuracy margin: %lf\n", epsilon);
    print_membership_matrix();
    printf("------------------------------------------------------------------------\n");
    return 0;
}