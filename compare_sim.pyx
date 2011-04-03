import numpy as np
cimport numpy as np

################################################################################
def compare_image_sims(img_list):
    ret = []

    for idxa in xrange(len(img_list)):
        for idxb in xrange(0, len(img_list)):
            sima, patha = img_list[idxa]
            simb, pathb = img_list[idxb]
            c = np.sum(np.absolute(np.subtract(sima, simb)))

            fp = (1.0 - (c / (255.0 * 1024.0 * 3.0)))

            if fp >= 0.98:
                ret.append( (fp, patha, pathb) )

    return ret


################################################################################
def compare_image_sims2(img_list):
    ret = []
    start = 0
    for idxa in xrange(len(img_list)):
        start += 1
        for idxb in xrange(start, len(img_list)):
            sima, patha = img_list[idxa]
            simb, pathb = img_list[idxb]
            c = np.sum(np.absolute(np.subtract(sima, simb)))

            fp = (1.0 - (c / (255.0 * 1024.0 * 3.0)))

            if fp >= 0.98:
                ret.append( (fp, patha, pathb) )

    return ret


################################################################################
def compare_image_sims3(img_list):
    cdef int start, idxa, idxb, length
    cdef double c, fp
#    cdef np.ndarray[np.float64_t, ndim=3] sima, simb

    ret = []
    start = 0
    length = len(img_list)

    for idxa in xrange(length):
        start += 1
        for idxb in xrange(start, length):
            sima, patha = img_list[idxa]
            simb, pathb = img_list[idxb]
            c = np.sum(np.absolute(np.subtract(sima, simb)))

            fp = (1.0 - (c / (255.0 * 1024.0 * 3.0)))

            if fp >= 0.98:
                ret.append( (fp, patha, pathb) )

    return ret


################################################################################
cpdef list compare_image_sims4(list img_list):
    cdef int start, idxa, idxb, length
    cdef double c, fp
#    cdef np.ndarray[np.float64_t, ndim=3] sima, simb
    cdef list ret

    DEF div = 255.0 * 1024.0 * 3.0

    ret = []
    start = 0
    length = len(img_list)

    for idxa in range(0, length):
        start += 1
        for idxb in range(start, length):
            sima, patha = img_list[idxa]
            simb, pathb = img_list[idxb]
            c = np.sum(np.absolute(np.subtract(sima, simb)))

            fp = (1.0 - (c / div))

            if fp >= 0.98:
                ret.append( (fp, patha, pathb) )

    return ret


################################################################################
cdef extern from "c_sim.h":
    void c_setup(int N)
    void c_add(int idx, double* a, char* b)

    void c_process1()
    void c_process2()

    void c_double_to_float()
    void c_process3()
    void c_process4()
    void c_teardown_floats()

    void c_teardown()

################################################################################
def setup(img_list):
    cdef int idx
    cdef np.ndarray[np.float64_t, ndim=3] sima

    c_setup(len(img_list))

    for idx in xrange(0, len(img_list)):
        sima, patha = img_list[idx]

        c_add( idx, <double*> sima.data, patha )

    c_double_to_float()
    c_process4()
    c_teardown_floats()

    c_teardown()
