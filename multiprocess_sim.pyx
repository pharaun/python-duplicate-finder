import numpy as np
cimport numpy as np
from multiprocessing import Pool

################################################################################
def sim_pool(lista):
    sima, patha, listb = lista
    ret = []

    for simb, pathb in listb:
        c = np.sum(np.absolute(np.subtract(sima, simb)))
        fp = (1.0 - (c / (255.0 * 1024.0 * 3.0)))

        if fp >= 0.98:
            ret.append( (fp, patha, pathb) )

    if not ret:
        return None
    else:
        return ret


################################################################################
cpdef list sim_pool2(tuple lista):
    cdef int idx, length
    cdef double c, fp
    cdef list ret, listb

    DEF div = 255.0 * 1024.0 * 3.0

    sima, patha, listb = lista
    ret = []
    length = len(listb)

    for idx in xrange(length):
        simb, pathb = listb[idx]
        c = np.sum(np.absolute(np.subtract(sima, simb)))
        fp = (1.0 - (c / div))

        if fp >= 0.98:
            ret.append( (fp, patha, pathb) )

    if not ret:
        return None
    else:
        return ret


################################################################################
def sim_pool_setup(img_list):
    ret = []
    listing = []

    i = 1
    for arr1, path1 in img_list:
        listing.append( (arr1, path1, img_list[i:]) )
        i += 1

    # Setup the pools?
    pool = Pool(processes=24)
    ret = pool.map(sim_pool2, listing)

    # Finish/close the pool
    pool.close()
    pool.join()

    # Flush out the None's
    ret = [item for item in ret if item != None]

    # Flatten the list
    ret = reduce(lambda x, y: x + y, ret)

    return ret


################################################################################
def sim_pool_setup2(list img_list):
    cdef int idx, length
    cdef list ret, listing

    ret = []
    listing = []
    length = len(img_list)

    for idx in xrange(length):
        arr, path = img_list[idx]
        listing.append( (arr, path, img_list[(idx+1):]))

    # Setup the pools?
    pool = Pool(processes=24)
    ret = pool.map(sim_pool2, listing)

    # Finish/close the pool
    pool.close()
    pool.join()

    # Flush out the None's
    ret = [item for item in ret if item != None]

    # Flatten the list
    ret = reduce(lambda x, y: x + y, ret)

    return ret


################################################################################
cdef list sim

def set_sim(n):
    global sim
    sim = n


################################################################################
def sim_pool3(lista):
    sima, patha, idx = lista
    length = len(sim)
    ret = []

    while idx < length:
        simb, pathb = sim[idx]
        c = np.sum(np.absolute(np.subtract(sima, simb)))
        fp = (1.0 - (c / (255.0 * 1024.0 * 3.0)))

        if fp >= 0.98:
            ret.append( (fp, patha, pathb) )

        idx += 1

    if not ret:
        return None
    else:
        return ret


################################################################################
cpdef list sim_pool4(tuple lista):
    cdef int length, idx
    cdef double c, fp
    cdef list ret

    DEF div = 255.0 * 1024.0 * 3.0

    sima, patha, idx = lista
    length = len(sim)
    ret = []

    while idx < length:
        simb, pathb = sim[idx]
        c = np.sum(np.absolute(np.subtract(sima, simb)))
        fp = (1.0 - (c / div))

        if fp >= 0.98:
            ret.append( (fp, patha, pathb) )

        idx += 1

    if not ret:
        return None
    else:
        return ret


################################################################################
def sim_pool_setup3():
    cdef int i, idx, length
    cdef list ret, listing

    ret = []
    listing = []
    length = len(sim)

    i = 1
    for idx in xrange(length):
        arr, path = sim[idx]
        listing.append( (arr, path, i) )
        i += 1

    # Setup the pools?
    pool = Pool(processes=24)
    ret = pool.map(sim_pool4, listing)

    # Finish/close the pool
    pool.close()
    pool.join()

    # Flush out the None's
    ret = [item for item in ret if item != None]

    # Flatten the list
    ret = reduce(lambda x, y: x + y, ret)

    return ret
