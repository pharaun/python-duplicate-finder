#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <time.h>

#include <omp.h>

#include <xmmintrin.h>
#include <pmmintrin.h>
#include <smmintrin.h>

// Each array == "3072 entry long - 3 * 32 * 32"
#define ARRAY_LENGTH 3072
// Similarity div == 255.0 * 1024.0 * 3.0 - Max value * element in array * rgb (3)
#define SIMILARITY_DIV 783360.0
#define SIMILARITY_THRESHOLD 0.98

// Number of ... element to create
#define COUNT 1000

// The arrays for each go here
const float** findex;
float* farray;

const uint8_t** uindex;
uint8_t* uarray;


//##############################################################################
// Create and destroy the test arrays
//##############################################################################
void create_arrays() {
    printf( "create_array\n");

    // The address of a block returned by malloc or realloc in the GNU system
    // is always a multiple of eight (or sixteen on 64-bit systems).
    findex = calloc(COUNT, sizeof(*findex));
    farray = calloc((COUNT * ARRAY_LENGTH), sizeof(*farray));

    uindex = calloc(COUNT, sizeof(*uindex));
    uarray = calloc((COUNT * ARRAY_LENGTH), sizeof(*uarray));

    if((findex == NULL) || (farray == NULL) || (uindex == NULL) || (uarray == NULL)) {
        printf("Out of memory\n");
    }

    // Init the random generator
    srand(time(NULL));

    int i, k;
    for(i = 0; i < COUNT; i++) {
        // File away this pointer as a index into the index array
        findex[i] = &farray[(i * ARRAY_LENGTH)];
        uindex[i] = &uarray[(i * ARRAY_LENGTH)];

        // Populate the farray with values
        for(k = 0; k < ARRAY_LENGTH; k++) {
	    // Generate random num here
	    uint8_t num = rand() % 256;

            farray[((i * ARRAY_LENGTH) + k)] = (float)num;
            uarray[((i * ARRAY_LENGTH) + k)] = num;
        }
    }

}

void destroy_arrays() {
    printf( "destroy_array\n");

    free(findex);
    free(farray);

    free(uindex);
    free(uarray);

    findex = NULL;
    farray = NULL;

    uindex = NULL;
    uarray = NULL;
}

//##############################################################################
// The float processing code
//##############################################################################
void c_float_process() {
    printf("c_float_process\n");

    int i, j, k;
    int dup = 0;

    #pragma omp parallel for shared(dup) private(j, i, k) schedule (dynamic)
    for(i = 0; i < COUNT; i++) {
        for(j = (i + 1); j < COUNT; j++) {

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
                }
            }
        }
    }

    printf("Duplicate found: %d\n", dup);
}

void sse_float_process() {
    printf("sse_float_process\n");

    int i, j, k;
    int dup = 0;

    #pragma omp parallel for shared(dup) private(j, i, k) schedule (dynamic)
    for(i = 0; i < COUNT; i++) {
        for(j = (i + 1); j < COUNT; j++) {

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
                #pragma omp critical
                {
                    dup += 1;
                }
            }
        }
    }

    printf("Duplicate found: %d\n", dup);
}

void sse_float_process2() {
    printf("sse_float_process2\n");

    int i, j, k;
    int dup = 0;

    #pragma omp parallel for shared(dup) private(j, i, k) schedule (dynamic)
    for(i = 0; i < COUNT; i++) {
        for(j = (i + 1); j < COUNT; j++) {

	    const __m128* sima = (__m128*) findex[i];
	    const __m128* simb = (__m128*) findex[j];

	    // Sign mask for abs - -0.0f = 1 << 31
	    const __m128 sign_mask = _mm_set1_ps(-0.0f);

	    // Init sum
	    __m128 sum1 = _mm_setzero_ps();
	    __m128 sum2 = _mm_setzero_ps();

            for(k = 0; k < ARRAY_LENGTH/4; k += 2) {
		sum1 += _mm_andnot_ps(sign_mask, _mm_sub_ps(sima[k], simb[k]));
		sum2 += _mm_andnot_ps(sign_mask, _mm_sub_ps(sima[k+1], simb[k+1]));
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
                }
            }
        }
    }

    printf("Duplicate found: %d\n", dup);
}

//##############################################################################
// The uint8 processing code
//##############################################################################
void c_uint8_process() {
    printf("c_uint8_process\n");

    int i, j, k;
    int dup = 0;

    #pragma omp parallel for shared(dup) private(j, i, k) schedule (dynamic)
    for(i = 0; i < COUNT; i++) {
        for(j = (i + 1); j < COUNT; j++) {

            const uint8_t* sima = uindex[i];
            const uint8_t* simb = uindex[j];

            int sum = 0;

            for(k = 0; k < ARRAY_LENGTH; k++) {
                sum += abs(sima[k] - simb[k]);
            }

            float fp = (1.0 - ((float)sum / SIMILARITY_DIV));

            if(fp >= SIMILARITY_THRESHOLD) {
                #pragma omp critical
                {
                    dup += 1;
                }
            }
        }
    }

    printf("Duplicate found: %d\n", dup);
}

void sse_uint8_process() {
    printf("sse_uint8_process\n");

    int i, j, k;
    int dup = 0;

    #pragma omp parallel for shared(dup) private(j, i, k) schedule (dynamic)
    for(i = 0; i < COUNT; i++) {
        for(j = (i + 1); j < COUNT; j++) {

            const uint8_t* sima = uindex[i];
            const uint8_t* simb = uindex[j];

	    // Init partial sums
	    __m128i vsum = _mm_setzero_si128();

	    for(k = 0; k < ARRAY_LENGTH; k += 16) {
		// Load 16 uint8_t from sima, simb
		__m128i va = _mm_load_si128((const __m128i*)&sima[k]);
		__m128i vb = _mm_load_si128((const __m128i*)&simb[k]);

		// Calc Sum Absolute Difference over the 16x uint8_t
		// 0, 0, 0, uint16 | 0, 0, 0, uint16
		__m128i vabsdiff = _mm_sad_epu8(va, vb);

		// Shuffle the high uint16 to the lower 64
		// 0, 0, 0, 0 | 0, uint16, 0, uint16
		vabsdiff = _mm_shuffle_epi32(vabsdiff, _MM_SHUFFLE(0, 2, 1, 3));
		// 0, 0, 0, 0 | 0, 0, uint16, uint16
		vabsdiff = _mm_shufflelo_epi16(vabsdiff, _MM_SHUFFLE(0, 2, 1, 3));

		// Convert the lower 4 uint16 to int32
		// 0, 0 | int32, int32
		vabsdiff = _mm_cvtepu16_epi32(vabsdiff);

		// Accumulate the two int32
		vsum = _mm_add_epi32(vsum, vabsdiff);
	    }

	    // Accumulate the partial sums into one
	    // 0, 0 | 0, int32
	    vsum = _mm_hadd_epi32(vsum, _mm_setzero_si128());

	    // Convert the signed int32 to float
	    __m128 fvsum = _mm_cvtepi32_ps(vsum);

	    // calc fvsum = fvsum / vdiv
	    fvsum = _mm_div_ps(fvsum, _mm_set1_ps(SIMILARITY_DIV));

	    // calc fvsum = 1.0 - fvsum
	    fvsum = _mm_sub_ps(_mm_set1_ps(1.0), fvsum);

	    // unload fvsum -> fp
	    float fp[4];
	    _mm_store_ps(&fp[0], fvsum);

	    if(fp[0] >= SIMILARITY_THRESHOLD) {
                #pragma omp critical
                {
                    dup += 1;
                }
            }
        }
    }

    printf("Duplicate found: %d\n", dup);
}

void sse_uint8_process2() {
    printf("sse_uint8_process2\n");

    int i, j, k;
    int dup = 0;

    #pragma omp parallel for shared(dup) private(j, i, k) schedule (dynamic)
    for(i = 0; i < COUNT; i++) {
        for(j = (i + 1); j < COUNT; j++) {

	    const __m128i* sima = (__m128i*) uindex[i];
	    const __m128i* simb = (__m128i*) uindex[j];

	    // Init partial sums
	    __m128i vsum = _mm_setzero_si128();

	    for(k = 0; k < ARRAY_LENGTH/16; k += 1) {
		vsum = _mm_add_epi32(vsum, _mm_sad_epu8(sima[k], simb[k]));
	    }

	    // Accumulate the partial sums into one
	    // 0, 0 | 0, int32
	    vsum = _mm_hadd_epi32(vsum, _mm_setzero_si128());
	    vsum = _mm_hadd_epi32(vsum, _mm_setzero_si128());

	    // Convert the signed int32 to float
	    __m128 fvsum = _mm_cvtepi32_ps(vsum);

	    // calc fvsum = fvsum / vdiv
	    fvsum = _mm_div_ps(fvsum, _mm_set1_ps(SIMILARITY_DIV));

	    // calc fvsum = 1.0 - fvsum
	    fvsum = _mm_sub_ps(_mm_set1_ps(1.0), fvsum);

	    // unload fvsum -> fp
	    float fp[4];
	    _mm_store_ps(&fp[0], fvsum);

	    if(fp[0] >= SIMILARITY_THRESHOLD) {
                #pragma omp critical
                {
                    dup += 1;
                }
            }
        }
    }

    printf("Duplicate found: %d\n", dup);
}

void sse_uint8_process3() {
    printf("sse_uint8_process3\n");

    int i, j, k;
    int dup = 0;

    #pragma omp parallel for shared(dup) private(j, i, k) schedule (dynamic)
    for(i = 0; i < COUNT; i++) {
        for(j = (i + 1); j < COUNT; j++) {

	    const __m128i* sima = (__m128i*) uindex[i];
	    const __m128i* simb = (__m128i*) uindex[j];

	    // Init sum
	    __m128i sum1 = _mm_setzero_si128();
	    __m128i sum2 = _mm_setzero_si128();

	    for(k = 0; k < ARRAY_LENGTH/16; k += 2) {
		sum1 = _mm_add_epi32(sum1, _mm_sad_epu8(sima[k], simb[k]));
		sum2 = _mm_add_epi32(sum2, _mm_sad_epu8(sima[k+1], simb[k+1]));
	    }

	    __m128i sum = _mm_add_epi32(sum1, sum2);

	    // Accumulate the partial sums into one
	    // 0, 0 | 0, int32
	    sum = _mm_hadd_epi32(sum, _mm_setzero_si128());
	    sum = _mm_hadd_epi32(sum, _mm_setzero_si128());

	    // Convert the signed int32 to float
	    __m128 fsum = _mm_cvtepi32_ps(sum);

	    // calc fvsum = fvsum / vdiv
	    fsum = _mm_div_ps(fsum, _mm_set1_ps(SIMILARITY_DIV));

	    // calc fvsum = 1.0 - fvsum
	    fsum = _mm_sub_ps(_mm_set1_ps(1.0), fsum);

	    // unload fsum -> fp
	    float fp[4];
	    _mm_store_ps(&fp[0], fsum);

	    if(fp[0] >= SIMILARITY_THRESHOLD) {
                #pragma omp critical
                {
                    dup += 1;
                }
            }
        }
    }

    printf("Duplicate found: %d\n", dup);
}


int main(int argc, const char** argv) {
    printf("Start\n");
    create_arrays();

    printf("Count: %d\n\n", COUNT);

    clock_t start, diff;
    int startepoch, diffepoch, msec;

    /*
    startepoch = time(NULL);
    start = clock();
    c_float_process();
    diff = clock() - start;
    diffepoch = time(NULL) - startepoch;
    msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Cpu time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);
    printf("User time taken %d seconds\n\n", diffepoch);

    startepoch = time(NULL);
    start = clock();
    sse_float_process();
    diff = clock() - start;
    diffepoch = time(NULL) - startepoch;
    msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Cpu time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);
    printf("User time taken %d seconds\n\n", diffepoch);

    startepoch = time(NULL);
    start = clock();
    sse_float_process2();
    diff = clock() - start;
    diffepoch = time(NULL) - startepoch;
    msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Cpu time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);
    printf("User time taken %d seconds\n\n", diffepoch);

    startepoch = time(NULL);
    start = clock();
    c_uint8_process();
    diff = clock() - start;
    diffepoch = time(NULL) - startepoch;
    msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Cpu time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);
    printf("User time taken %d seconds\n\n", diffepoch);

    startepoch = time(NULL);
    start = clock();
    sse_uint8_process();
    diff = clock() - start;
    diffepoch = time(NULL) - startepoch;
    msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Cpu time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);
    printf("User time taken %d seconds\n\n", diffepoch);
    */

    startepoch = time(NULL);
    start = clock();
    sse_uint8_process2();
    diff = clock() - start;
    diffepoch = time(NULL) - startepoch;
    msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Cpu time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);
    printf("User time taken %d seconds\n\n", diffepoch);

    startepoch = time(NULL);
    start = clock();
    sse_uint8_process3();
    diff = clock() - start;
    diffepoch = time(NULL) - startepoch;
    msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Cpu time taken %d seconds %d milliseconds\n", msec/1000, msec%1000);
    printf("User time taken %d seconds\n\n", diffepoch);

    destroy_arrays();
    printf("End\n");
}
