#define SSE

#include "c_sim.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <omp.h>

#include <xmmintrin.h>
#include <pmmintrin.h>
#include <smmintrin.h>

// Each array == "3072 entry long - 3 * 32 * 32"
#define ARRAY_LENGTH 3072
// Similarity div == 255.0 * 1024.0 * 3.0 - Max value * element in array * rgb (3)
#define SIMILARITY_DIV 783360.0
#define SIMILARITY_THRESHOLD 0.98

// Global vars
//append_to_list()
struct node* head = NULL;
struct node* last = NULL;
int list_length = 0;

//c_get_similarity_next()
struct node* cur = NULL;
int done = 0;

// Rest of the app
const double** array;
const char** path;
int length;

// Floats part of the app
const float** findex;
float* farray;


//##############################################################################
// Simple Append only linked list
//##############################################################################
int append_to_list(double fp, int patha, int pathb) {
    // allocate the new struct
    struct node* newnode = malloc(sizeof(struct node));

    if(newnode == NULL) {
	printf("Out of memory\n");
    }

    // set the data
    newnode->fp = fp;
    newnode->patha = patha;
    newnode->pathb = pathb;

    newnode->next = NULL;

    // Check to see if list is empty
    if(head == NULL) {
	head = newnode;
	last = newnode;

	list_length += 1;
	return 0;
    }

    // List is not empty
    last->next = newnode;
    last = last->next;

    list_length += 1;

    return 0;
}


int free_list() {
    if(head == NULL) {
	return 0;
    }

    struct node* prev = head;
    struct node* next = head->next;

    int count = 0;
    while(prev != NULL) {
	free(prev);
	count += 1;

	prev = next;
	if(prev != NULL) {
	    next = prev->next;
	}
    }

    printf("freed: %d\n", count);

    head = NULL;
    last = NULL;
    list_length = 0;

    return 0;
}

//##############################################################################
// Return the values out of the linked list
//##############################################################################
int c_get_similarity_next(struct similar* a) {
    if(head == NULL) {
	return 0;
    } else if(done) {
	return 0;
    } else if(cur == NULL) {
	cur = head;
    }

    // Fill the struct
    a->fp = cur->fp;
    a->patha = path[cur->patha];
    a->pathb = path[cur->pathb];

    struct node* next = cur->next;
    if(next != NULL) {
	cur = next;
    } else {
	done = 1;
    }

    return 1;
}

//##############################################################################
// Setup and cleanup code
//##############################################################################
void c_setup( int N ) {
    printf( "c_setup: %d\n", N);

    length = N;
    array = calloc(N, sizeof(*array));
    path = calloc(N, sizeof(*path));

    if((array == NULL) || (path == NULL)) {
	printf("Out of memory\n");
    }
}

void c_add( int idx, const double a[], const char* b ) {
    if((idx >= 0) && (idx < length)) {
	array[idx] = &a[0];
	path[idx] = b;
    }
}

void c_teardown() {
    printf( "c_teardown\n" );

    free(array);
    free(path);

    array = NULL;
    path = NULL;

    free_list();
}

//##############################################################################
// Double -> Floats here
//##############################################################################
void c_double_to_float() {
    printf( "c_double_to_float:\n");

    // The address of a block returned by malloc or realloc in the GNU system
    // is always a multiple of eight (or sixteen on 64-bit systems). 
    findex = calloc(length, sizeof(*findex));
    farray = calloc((length * ARRAY_LENGTH), sizeof(*farray));

    if((findex == NULL) || (farray == NULL)) {
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

//##############################################################################
// Inlined abs
//##############################################################################
inline __m128 abs_ps(__m128 x) {
    const __m128 sign_mask = _mm_set1_ps(-0.0f); // -0.0f = 1 << 31
    return _mm_andnot_ps(sign_mask, x);
}

inline __m128d abs_pd(__m128d x) {
    const __m128d sign_mask = _mm_set1_pd(-0.0); // -0.0 = 1 << 63
    return _mm_andnot_pd(sign_mask, x); // !sign_mask & x
}

//##############################################################################
// Setup and cleanup code
//##############################################################################
void c_process() {
    printf( "c_process\n");

    int i, j, k;
    int dup = 0;

#ifdef DOUBLE
#ifndef SSE
    #pragma omp parallel for shared(dup) private(j, i, k) schedule (dynamic)
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
		#pragma omp critical
		{
		    dup += 1;
		    append_to_list(fp, i, j);
		}
	    }
	}
    }
#else
    #pragma omp parallel for shared(dup) private(j, i, k) schedule (dynamic)
    for(i = 0; i < length; i++) {
	for(j = (i + 1); j < length; j++) {

	    const __m128d* sima = (__m128d*) array[i];
	    const __m128d* simb = (__m128d*) array[j];

	    // Init sums
	    __m128d sum1 = _mm_setzero_pd();
	    __m128d sum2 = _mm_setzero_pd();

	    for(k = 0; k < ARRAY_LENGTH/2; k += 2) {
		sum1 += abs_pd(_mm_sub_pd(sima[k], simb[k]));
		sum2 += abs_pd(_mm_sub_pd(sima[k+1], simb[k+1]));
	    }

	    __m128d sum = _mm_add_pd(sum1, sum2);

	    // Accumulate the partial sums into one
	    sum = _mm_hadd_pd(sum, _mm_setzero_pd());

	    // calc vsum = vsum / vdiv
	    sum = _mm_div_sd(sum, _mm_set1_pd(SIMILARITY_DIV));

	    // calc vsum = 1.0 - vsum
	    sum = _mm_sub_pd(_mm_set1_pd(1.0), sum);

	    // Unload vsum -> fp
	    double fp[2];
	    _mm_store_sd(&fp[0], sum);

	    if(fp[0] >= SIMILARITY_THRESHOLD) {
		#pragma omp critical
		{
		    dup += 1;
		    append_to_list(fp[0], i, j);
		}
	    }
	}
    }
#endif
#else
    // Convert the doubles to float
    c_double_to_float();
#ifndef SSE
    #pragma omp parallel for shared(dup) private(j, i, k) schedule (dynamic)
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
		#pragma omp critical
		{
		    dup += 1;
		    append_to_list((double)fp, i, j);
		}
	    }
	}
    }
#else
    #pragma omp parallel for shared(dup) private(j, i, k) schedule (dynamic)
    for(i = 0; i < length; i++) {
	for(j = (i + 1); j < length; j++) {

	    const __m128* sima = (__m128*) findex[i];
	    const __m128* simb = (__m128*) findex[j];

	    // Init sums
	    __m128 sum1 = _mm_setzero_ps();
	    __m128 sum2 = _mm_setzero_ps();

	    for(k = 0; k < ARRAY_LENGTH/4; k += 2) {
		sum1 += abs_ps(_mm_sub_ps(sima[k], simb[k]));
		sum2 += abs_ps(_mm_sub_ps(sima[k+1], simb[k+1]));
	    }

	    __m128 sum = _mm_add_ps(sum1, sum2);

	    // Accumulate the partial sums into one
	    sum = _mm_hadd_ps(sum, _mm_setzero_ps());
	    sum = _mm_hadd_ps(sum, _mm_setzero_ps());

	    // calc vsum = vsum / vdiv
	    sum = _mm_div_ps(sum, _mm_set1_ps(SIMILARITY_DIV));

	    // calc vsum = 1.0 - vsum
	    sum = _mm_sub_ps(_mm_set1_ps(1.0), sum);

	    // unload vsum -> fp
	    float fp[4];
	    _mm_store_ps(&fp[0], sum);

	    if(fp[0] >= SIMILARITY_THRESHOLD) {
		#pragma omp critical
		{
		    dup += 1;
		    append_to_list((double)fp[0], i, j);
		}
	    }
	}
    }
#endif

    // Clean up
    c_teardown_floats();
#endif
    printf("dup: %d - list-length: %d\n", dup, list_length );
}
