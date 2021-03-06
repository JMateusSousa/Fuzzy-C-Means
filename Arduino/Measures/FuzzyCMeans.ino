/* * @(#)FuzzyCMeans.c
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


#include <avr/pgmspace.h>

#define NUM_DATA_POINTS 60
#define NUM_CLUSTERS 3
#define NUM_DIMENSIONS 13

float epsilon = 0.0001;
float fuzziness = 2.00;
float degree_of_memb[NUM_DATA_POINTS][NUM_CLUSTERS];
float cluster_centre[NUM_CLUSTERS][NUM_DIMENSIONS];
float t[NUM_DATA_POINTS][NUM_CLUSTERS];
unsigned long start_time, end_time;

static const float data_point[][NUM_DIMENSIONS] PROGMEM = {
{14.23,1.71,2.43,15.6,127,2.8,3.06,0.28,2.29,5.64,1.04,3.92,1065},
{13.2,1.78,2.14,11.2,100,2.65,2.76,0.26,1.28,4.38,1.05,3.4,1050},
{13.16,2.36,2.67,18.6,101,2.8,3.24,0.3,2.81,5.68,1.03,3.17,1185},
{14.37,1.95,2.5,16.8,113,3.85,3.49,0.24,2.18,7.8,0.86,3.45,1480},
{13.24,2.59,2.87,21,118,2.8,2.69,0.39,1.82,4.32,01.04,2.93,735},
{14.2,1.76,2.45,15.2,112,3.27,3.39,0.34,1.97,6.75,01.05,2.85,1450},
{14.39,1.87,2.45,14.6,96,2.5,2.52,0.3,1.98,5.25,01.02,3.58,1290},
{14.06,2.15,2.61,17.6,121,2.6,2.51,0.31,1.25,05.05,01.06,3.58,1295},
{14.83,1.64,2.17,14,97,2.8,2.98,0.29,1.98,5.2,01.08,2.85,1045},
{13.86,1.35,2.27,16,98,2.98,3.15,0.22,1.85,7.22,01.01,3.55,1045},
{14.1,2.16,2.3,18,105,2.95,3.32,0.22,2.38,5.75,1.25,3.17,1510},
{14.12,1.48,2.32,16.8,95,2.2,2.43,0.26,1.57,5,1.17,2.82,1280},
{13.75,1.73,2.41,16,89,2.6,2.76,0.29,1.81,5.6,1.15,2.9,1320},
{14.75,1.73,2.39,11.4,91,3.1,3.69,0.43,2.81,5.4,1.25,2.73,1150},
{14.38,1.87,2.38,12,102,3.3,3.64,0.29,2.96,7.5,1.2,3,1547},
{13.63,1.81,2.7,17.2,112,2.85,2.91,0.3,1.46,7.3,1.28,2.88,1310},
{14.3,1.92,2.72,20,120,2.8,3.14,0.33,1.97,6.2,01.07,2.65,1280},
{13.83,1.57,2.62,20,115,2.95,3.4,0.4,1.72,6.6,1.13,2.57,1130},
{14.19,1.59,2.48,16.5,108,3.3,3.93,0.32,1.86,8.7,1.23,2.82,1680},
{13.64,3.1,2.56,15.2,116,2.7,03.03,0.17,1.66,5.1,0.96,3.36,845},
{12.21,1.19,1.75,16.8,151,1.85,1.28,0.14,2.5,2.85,1.28,03.07,718},
{12.29,1.61,2.21,20.4,103,1.1,01.02,0.37,1.46,03.05,906,1.82,870},
{13.86,1.51,2.67,25,86,2.95,2.86,0.21,1.87,3.38,1.36,3.16,410},
{13.49,1.66,2.24,24,87,1.88,1.84,0.27,01.03,3.74,0.98,2.78,472},
{12.99,1.67,2.6,30,139,3.3,2.89,0.21,1.96,3.35,1.31,3.5,985},
{11.96,01.09,2.3,21,101,3.38,2.14,0.13,1.65,3.21,0.99,3.13,886},
{11.66,1.88,1.92,16,97,1.61,1.57,0.34,1.15,3.8,1.23,2.14,428},
{13.03,0.9,1.71,16,86,1.95,02.03,0.24,1.46,4.6,1.19,2.48,392},
{11.84,2.89,2.23,18,112,1.72,1.32,0.43,0.95,2.65,0.96,2.52,500},
{12.7,3.87,2.4,23,101,2.83,2.55,0.43,1.95,2.57,1.19,3.13,463},
{12,0.92,2,19,86,2.42,2.26,0.3,1.43,2.5,1.38,3.12,278},
{12.72,1.81,2.2,18.8,86,2.2,2.53,0.26,1.77,3.9,1.16,3.14,714},
{12.08,1.13,2.51,24,78,2,1.58,0.4,1.4,2.2,1.31,2.72,630},
{13.05,3.86,2.32,22.5,85,1.65,1.59,0.61,1.62,4.8,0.84,02.01,515},
{11.84,0.89,2.58,18,94,2.2,2.21,0.22,2.35,03.05,0.79,03.08,520},
{12.67,0.98,2.24,18,99,2.2,1.94,0.3,1.46,2.62,1.23,3.16,450},
{12.16,1.61,2.31,22.8,90,1.78,1.69,0.43,1.56,2.45,1.33,2.26,495},
{11.65,1.67,2.62,26,88,1.92,1.61,0.4,1.34,2.6,1.36,3.21,562},
{11.64,02.06,2.46,21.6,84,1.95,1.69,0.48,1.35,2.8,1,2.75,680},
{12.08,1.33,2.3,23.6,70,2.2,1.59,0.42,1.38,1.74,01.07,3.21,625},
{12.45,3.03,2.64,27,97,1.9,0.58,0.63,1.14,7.5,0.67,1.73,880},
{14.34,1.68,2.7,25,98,2.8,1.31,0.53,2.7,13,0.57,1.96,660},
{13.48,1.67,2.64,22.5,89,2.6,1.1,0.52,2.29,11.75,0.57,1.78,620},
{12.36,3.83,2.38,21,88,2.3,0.92,0.5,01.04,7.65,0.56,1.58,520},
{13.69,3.26,2.54,20,107,1.83,0.56,0.5,0.8,5.88,0.96,1.82,680},
{12.85,3.27,2.58,22,106,1.65,0.6,0.6,0.96,5.58,0.87,2.11,570},
{12.96,3.45,2.35,18.5,106,1.39,0.7,0.4,0.94,5.28,0.68,1.75,675},
{13.78,2.76,2.3,22,90,1.35,0.68,0.41,01.03,9.58,0.7,1.68,615},
{13.73,4.36,2.26,22.5,88,1.28,0.47,0.52,1.15,6.62,0.78,1.75,520},
{13.45,3.7,2.6,23,111,1.7,0.92,0.43,1.46,10.68,0.85,1.56,695},
{12.82,3.37,2.3,19.5,88,1.48,0.66,0.4,0.97,10.26,0.72,1.75,685},
{13.58,2.58,2.69,24.5,105,1.55,0.84,0.39,1.54,8.66,0.74,1.8,750},
{13.4,4.6,2.86,25,112,1.98,0.96,0.27,1.11,8.5,0.67,1.92,630},
{12.2,3.03,2.32,19,96,1.25,0.49,0.4,0.73,5.5,0.66,1.83,510},
{12.77,2.39,2.28,19.5,86,1.39,0.51,0.48,0.64,9.899999,0.57,1.63,470},
{14.16,2.51,2.48,20,91,1.68,0.7,0.44,1.24,9.7,0.62,1.71,660},
{13.71,5.65,2.45,20.5,95,1.68,0.61,0.52,01.06,7.7,0.64,1.74,740},
{13.4,3.91,2.48,23,102,1.8,0.75,0.43,1.41,7.3,0.7,1.56,750},
{13.27,4.28,2.26,20,120,1.59,0.69,0.43,1.35,10.2,0.59,1.56,835},
{13.17,2.59,2.37,20,120,1.65,0.68,0.53,1.46,9.3,0.6,1.62,840},
{14.13,4.1,2.74,24.5,96,2.05,0.76,0.56,1.35,9.2,0.61,1.6,560}};

void setup() {
  Serial.begin(9600);
  fuzzy_c_means();
}

void loop() { 
}

void init_matrix() {
  int i, j, r, rval, aux;
    float s;
  for (i = 0; i < NUM_DATA_POINTS; i++) {
        s = 0.0;
        r = 100;
    aux = i;
        for (j = 1; j < NUM_CLUSTERS; j++) {
            rval =  rand() % (r + 1);
            r -= rval;
            degree_of_memb[i][j] = rval / 100.0;
            s += degree_of_memb[i][j];
      aux++;
        }
        degree_of_memb[i][0] = 1.0 - s;
    }
}

int calculate_centre_vectors() {
    int i, j, k;
    float numerator, denominator;
    for (i = 0; i < NUM_DATA_POINTS; i++) {
        for (j = 0; j < NUM_CLUSTERS; j++) {
            t[i][j] = pow(degree_of_memb[i][j], fuzziness);
        }
    }
    for (j = 0; j < NUM_CLUSTERS; j++) {
        for (k = 0; k < NUM_DIMENSIONS; k++) {
            numerator = 0.0;
            denominator = 0.0;
            for (i = 0; i < NUM_DATA_POINTS; i++) {
                numerator += t[i][j] * pgm_read_float_near(&data_point[i][k]);
                denominator += t[i][j];
            }
            cluster_centre[j][k] = numerator / denominator;
        }
    }
    return 0;
}

float get_norm(int i, int j) {
    int k;
    float sum = 0.0;
    for (k = 0; k < NUM_CLUSTERS; k++) {
        sum += pow(pgm_read_float_near(&data_point[i][k]) - cluster_centre[j][k], 2);
    }
    return sqrt(sum);
}

float get_new_value(int i, int j) {
    int k;
    float t, p, sum;
    sum = 0.0;
    p = 2 / (fuzziness - 1);
    for (k = 0; k < NUM_CLUSTERS; k++) {
        t = get_norm(i, j) / get_norm(i, k);
        t = pow(t, p);
        sum += t;
    }
    return 1.0 / sum;
}

float update_degree_of_membership() {
    int i, j;
    float new_uij;
    float max_diff = 0.0, diff;
    for (j = 0; j < NUM_CLUSTERS; j++) {
        for (i = 0; i < NUM_DATA_POINTS; i++) {
            new_uij = get_new_value(i, j);
            diff = new_uij - degree_of_memb[i][j];
            if (diff > max_diff)
                max_diff = diff;
            degree_of_memb[i][j] = new_uij;
        }
    }
    return max_diff;
}

void fuzzy_c_means() {
    float max_diff;
    start_time = micros();
  init_matrix();
  int i, j; 
    int clusters[3] = {0};
    do {
        calculate_centre_vectors();
        max_diff = update_degree_of_membership();
    } while (max_diff > epsilon);
    end_time = micros();
  Serial.println(start_time - end_time);
  for (i = 0; i < NUM_DATA_POINTS; i++) {
        float bigger = 0;
        char index = 0;
        for (j = 0; j < NUM_CLUSTERS; j++) {
            if (bigger < degree_of_memb[i][j]) {
                bigger = degree_of_memb[i][j];
                index = j;
            }
        }
        clusters[index]++;
  }
  Serial.print("grupo 01: ");
  Serial.println(clusters[0]);
  Serial.print("grupo 02: ");
  Serial.println(clusters[1]);
  Serial.print("grupo 03: ");
  Serial.println(clusters[2]);
}
