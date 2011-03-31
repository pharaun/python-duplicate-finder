#include "c_sim.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <omp.h>

const double** array;
const char** path;
int length;

void c_setup( int N ) {
    printf( "c_setup: %d\n", N);

    length = N;
    array = malloc(N * sizeof(*array));
    path = malloc(N * sizeof(*path));

    if ((array == NULL) || (path == NULL)) {
	printf("Out of memory\n");
    }
}

void c_add( int idx, const double a[], const char* b ) {
    if ((idx >= 0) && (idx < length)) {
	array[idx] = &a[0];
	path[idx] = b;
    }
}

void c_process() {
    // Each array == "3072 entry long - 3 * 32 * 32"
    printf( "c_process\n");

    int i, j, k;
    int dup = 0;

    #pragma omp parallel for shared(dup, array, path) private(j, i, k)
    for(i = 0; i < length; i++) {
	for(j = (i + 1); j < length; j++) {
	    
	    const double* sima = array[i];
	    const double* simb = array[j];

	    double sum = 0.0;

	    for(k = 0; k < 3072; k++) {
		sum += fabs(sima[k] - simb[k]);
	    }

	    double fp = (1.0 - (sum / (255.0 * 1024.0 * 3.0)));

	    if(fp >= 0.98) {
//		printf("%f - %s - %s\n", fp, path[i], path[j]);
		#pragma omp atomic
		dup += 1;
	    }
	}
    }

    printf("dup: %d\n", dup);
}

void c_teardown() {
    printf( "c_teardown\n" );

    free(array);
    free(path);

    array = NULL;
    path = NULL;
}
