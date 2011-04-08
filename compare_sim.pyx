import numpy as np
cimport numpy as np

################################################################################
cdef extern from "c_sim.h":
    cdef struct similar:
        double fp
        char* patha
        char* pathb

    int c_get_similarity_next(similar* a)

    void c_setup(int N)
    void c_add(int idx, double* a, char* b)
    void c_process()
    void c_teardown()

################################################################################
def compare(img_list):
    cdef int idx, length
    cdef np.ndarray[np.float64_t, ndim=3] sima
    cdef list ret
    cdef similar tmp

    length = len(img_list)

    # Setup the indexing data structure & etc
    c_setup(length)

    # Transfer it to the C module
    for idx in xrange(0, length):
        sima, patha = img_list[idx]
        c_add( idx, <double*> sima.data, patha )

    # Process it
    c_process()

    # Fetch the information
    ret = []

    while (c_get_similarity_next(&tmp) != 0):
        ret.append( (tmp.patha, tmp.pathb, tmp.fp) )

    # Clean it all up
    c_teardown()

    return ret
