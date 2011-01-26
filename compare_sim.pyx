import numpy as np

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
