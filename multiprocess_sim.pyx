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
def sim_pool_setup(img_list):
    ret = []
    listing = []

    i = 1
    for arr1, path1 in img_list:
        listing.append( (arr1, path1, img_list[i:]) )
        i += 1

    # Setup the pools?
    pool = Pool(processes=24)
    ret = pool.map(sim_pool, listing)

    # Finish/close the pool
    pool.close()
    pool.join()

    # Flush out the None's
    ret = [item for item in ret if item != None]

    # Flatten the list
    ret = reduce(lambda x, y: x + y, ret)

    return ret
