#include "c_sim.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <omp.h>

#include <xmmintrin.h>
#include <pmmintrin.h>

// Each array == "3072 entry long - 3 * 32 * 32"
#define ARRAY_LENGTH 3072
// Similarity div == 255.0 * 1024.0 * 3.0 - Max value * element in array * rgb (3)
#define SIMILARITY_DIV 783360.0
#define SIMILARITY_THRESHOLD 0.98

// Macros to hint on branching
#define likely(x)       __builtin_expect((x),1)
#define unlikely(x)     __builtin_expect((x),0)

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
    printf( "c_process1\n");

    int i, j, k;
    int dup = 0;

   #pragma omp parallel for shared(dup, array, path) private(j, i, k) schedule (dynamic)
    for(i = 0; i < length; i++) {
	for(j = (i + 1); j < length; j++) {
	    
	    const double* sima = array[i];
	    const double* simb = array[j];

	    double sum = 0.0;

	    for(k = 0; k < ARRAY_LENGTH; k++) {
		sum += fabs(sima[k] - simb[k]);
	    }

	    double fp = (1.0 - (sum / SIMILARITY_DIV));

	    if(fp >= SIMILARITY_THRESHOLD) {
//		printf("%f - %s - %s\n", fp, path[i], path[j]);
		#pragma omp atomic
		dup += 1;
	    }
	}
    }

    printf("dup: %d\n", dup);
}

void c_process2() {
    printf( "c_process2\n");

    int i, j, k;
    int dup = 0;

    #pragma omp parallel for shared(dup, array, path) private(j, i, k) schedule (dynamic)
    for(i = 0; i < length; i++) {
	for(j = (i + 1); j < length; j++) {

	    const double* sima = array[i];
	    const double* simb = array[j];

	    // Init partial sums
	    __m128d vsum = _mm_set1_pd(0.0);

	    for(k = 0; k < ARRAY_LENGTH; k += 2) {
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
	    vsum = _mm_div_sd(vsum, _mm_set1_pd(SIMILARITY_DIV));

	    // calc vsum = 1.0 - vsum
	    vsum = _mm_sub_pd(_mm_set1_pd(1.0), vsum);

	    // Unload vsum -> fp
	    double fp[2];
	    _mm_store_sd(&fp[0], vsum);

	    if(fp[0] >= SIMILARITY_THRESHOLD) {
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


//##############################################################################
// Double -> Floats here
//##############################################################################

const float** findex;
float* farray;

void c_double_to_float() {
    printf( "c_double_to_float:\n");
    
    findex = calloc(length, sizeof(*findex));
    farray = calloc((length * ARRAY_LENGTH), sizeof(*farray));

    if ((findex == NULL) || (farray == NULL)) {
	printf("Out of memory\n");
    }

    // Start the conversion here
    int i, k;

    for(i = 0; i < length; i++) {
	const double* sim = array[i];

	// File away this pointer as a index into the index array
	findex[i] = &farray[(i * ARRAY_LENGTH)];

	// Populate the farray with values
	for(k = 0; k < ARRAY_LENGTH; k++) {
	    farray[((i * ARRAY_LENGTH) + k)] = (float)sim[k];
	}
    }
}

void c_teardown_floats() {
    printf( "c_teardown_floats\n" );

    free(findex);
    free(farray);

    findex = NULL;
    farray = NULL;
}

void c_process3() {
    printf( "c_process3\n");

    int i, j, k;
    int dup = 0;

   #pragma omp parallel for shared(dup, array, path) private(j, i, k) schedule (dynamic)
    for(i = 0; i < length; i++) {
	for(j = (i + 1); j < length; j++) {
	    
	    const float* sima = findex[i];
	    const float* simb = findex[j];

	    float sum = 0.0;

	    for(k = 0; k < ARRAY_LENGTH; k++) {
		sum += fabs(sima[k] - simb[k]);
	    }

	    float fp = (1.0 - (sum / SIMILARITY_DIV));

	    if(fp >= SIMILARITY_THRESHOLD) {
//		printf("%f - %s - %s\n", fp, path[i], path[j]);
		#pragma omp atomic
		dup += 1;
	    }
	}
    }

    printf("dup: %d\n", dup);
}

void c_process4() {
    printf( "c_process4\n");

    int i, j, k;
    int dup = 0;

   #pragma omp parallel for shared(dup, array, path) private(j, i, k) schedule (dynamic)
    for(i = 0; i < length; i++) {
	for(j = (i + 1); j < length; j++) {
	    
	    const float* sima = findex[i];
	    const float* simb = findex[j];

	    // Init partial sums
	    __m128 vsum = _mm_setzero_ps();

	    for(k = 0; k < ARRAY_LENGTH; k += 4) {
		// Load 4 floats from sima, simb
		__m128 va = _mm_load_ps(&sima[k]);
		__m128 vb = _mm_load_ps(&simb[k]);

		// calc diff = sima - simb
		__m128 vdiff = _mm_sub_ps(va, vb);

		// calc neg diff = 0.0 - diff
		__m128 vnegdiff = _mm_sub_ps(_mm_setzero_ps(), vdiff);

		// calc abs diff = max(diff, - diff)
		__m128 vabsdiff = _mm_max_ps(vdiff, vnegdiff);

		// accumulate two partial sums
		vsum = _mm_add_ps(vsum, vabsdiff);
	    }

	    // Accumulate the partial sums into one
	    vsum = _mm_hadd_ps(vsum, _mm_setzero_ps());
	    vsum = _mm_hadd_ps(vsum, _mm_setzero_ps());

	    // calc vsum = vsum / vdiv
	    vsum = _mm_div_ps(vsum, _mm_set1_ps(SIMILARITY_DIV));

	    // calc vsum = 1.0 - vsum
	    vsum = _mm_sub_ps(_mm_set1_ps(1.0), vsum);

	    // unload vsum -> fp
	    float fp[4];
	    _mm_store_ps(&fp[0], vsum);

	    if(fp[0] >= SIMILARITY_THRESHOLD) {
//		printf("%f - %s - %s\n", fp, path[i], path[j]);
		#pragma omp atomic
		dup += 1;
	    }
	}
    }

    printf("dup: %d\n", dup);
}
