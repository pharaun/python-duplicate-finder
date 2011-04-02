#include "c_sim.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <omp.h>

#include <xmmintrin.h>
#include <pmmintrin.h>

const double** array;
const char** path;
int length;

void c_setup( int N ) {
    printf( "c_setup: %d\n", N);

    length = N;
//    array = malloc(N * sizeof(*array));
    array = calloc(N, sizeof(*array));
//    path = malloc(N * sizeof(*path));
    path = calloc(N, sizeof(*path));

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

void c_process1() {
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

void c_process2() {
    // Each array == "3072 entry long - 3 * 32 * 32"
    printf( "c_process\n");

    int i, j, k;
    int dup = 0;

    for(i = 0; i < length; i++) {
	for(j = (i + 1); j < length; j++) {

	    const double* sima = array[i];
	    const double* simb = array[j];

	    // Init partial sums
	    __m128d vsum = _mm_set1_pd(0.0);

	    for(k = 0; k < 3072; k += 2) {
		// Load 2 doubles from sima, simb
		__m128d va = _mm_load_pd(&sima[k]);
		__m128d vb = _mm_load_pd(&simb[k]);

		// Calc diff = sima - simb
		__m128d vdiff = _mm_sub_pd(va, vb);

		// calc neg diff = 0.0 - diff
		__m128d vnegdiff = _mm_sub_pd(_mm_set1_pd(0.0), vdiff);

		// calc abs diff = max(diff, - diff)
		__m128d vabsdiff = _mm_max_pd(vdiff, vnegdiff);

		// accumulate two partial sums
		vsum = _mm_add_pd(vsum, vabsdiff);
	    }

	    // Accumulate the partial sums into one
	    vsum = _mm_hadd_pd(vsum, _mm_set1_pd(0.0));

	    // calc vsum = vsum / vdiv
	    vsum = _mm_div_sd(vsum, _mm_set1_pd(255.0 * 1024.0 * 3.0));

	    // calc vsum = 1.0 - vsum
	    vsum = _mm_sub_pd(_mm_set1_pd(1.0), vsum);

	    // Unload vsum -> fp
	    double fp[2];
	    _mm_store_sd(&fp[0], vsum);

	    if(fp[0] >= 0.98) {
//		printf("%f - %s - %s\n", fp, path[i], path[j]);
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
