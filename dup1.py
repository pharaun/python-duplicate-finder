#!/usr/bin/env python
"""
Python Duplicate picture finder, this mainly tests
various stuff related to the duplicate finding, such as
file/directory access, imhdr(identification) and so forth
"""

import os
import sys
import imghdr
import hashlib
import Image
import numpy as np
from multiprocessing import Pool
from optparse import OptionParser
from guppy import hpy


################################################################################
def calc_image_stats(img_list):
    count = {}
    total = 0

    for root, filename in img_list:
        try:
            im = Image.open(os.path.join(root, filename))
            count[(im.mode, im.format)] = count.setdefault((im.mode, im.format), 0) + 1
            total += 1
        except IOError:
            # Ignore bad images
            1+1

    # Print results
    mode = ['1', 'P', 'L', 'LA', 'RGB', 'RGBA', 'RGBX', 'CMYK']
    imgf = ['BMP', 'PNG', 'GIF', 'JPEG', 'TIFF']

    for m in mode:
        print m
        for t in imgf:
            try:
                v = count[(m, t)]
                print "\t" + t + "\t - %.2f" % (100.0 * (float(v) / float(total))) + "%\t - " + str(v)
            except KeyError:
                # Ignore
                1+1

    print "---"
    print "TOTAL\t\t - 100.00%\t - " + str(total)


################################################################################
def generate_img_list(rootdir, exclude):
    img = []

    for root, dirs, files in os.walk(rootdir):
        if exclude in dirs:
            dirs.remove(exclude)

        for filename in files:
            path = os.path.abspath(os.path.join(root, filename))
            what = imghdr.what(path)

            if what != None:
                img.append( (root, filename) )
    return img


################################################################################
def generate_image_sim(path):
    mode = ('RGB',)

    try:
        with open(path) as f:
            image = Image.open(f)
            if not (image.mode in mode):
                return None

            width, height = image.size

            x_inc = width / 32
            y_inc = height / 32

            # For now don't deal with small images
            if x_inc < 1:
                return None
            elif y_inc < 1:
                return None

            image_pixel = np.asarray(image)
            img_sum = np.zeros( (32, 32, 3) )

            # Start looping through it
            j = 0
            for ys in range(32):
                i = 0
                for xs in range(32):
                    # Sanity check
                    try:
                        r = np.mean(image_pixel[i:i+y_inc, j:j+x_inc, 0:1])
                        g = np.mean(image_pixel[i:i+y_inc, j:j+x_inc, 1:2])
                        b = np.mean(image_pixel[i:i+y_inc, j:j+x_inc, 2:3])
                    except IndexError:
                        return None

                    img_sum[ys:, xs:, 0:1] = r
                    img_sum[ys:, xs:, 1:2] = g
                    img_sum[ys:, xs:, 2:3] = b

                    i += x_inc
                j += y_inc

            return (np.rint(np.nan_to_num(img_sum)), path)

    except IOError:
        return None


################################################################################
def generate_sim_data(img_list):
    img = []

    filepath = []
    for root, filename in img_list:
        filepath.append( os.path.join(root, filename))

    # Setup the pools?
    pool = Pool(processes=24)
    img = pool.map(generate_image_sim, filepath)

    # Finish/close the pool
    pool.close()
    pool.join()

    # Flush out the None's
    img = [item for item in img if item != None]
    return img


################################################################################
def image_sim_compare_pool(lista):
    simpatha, patha, listb = lista
    ret = []

    for simpathb, pathb in listb:
        sima = np.load(simpatha)
        simb = np.load(simpathb)
        c = np.sum(np.absolute(np.subtract(sima, simb)))
        fp = (1.0 - (c / (255.0 * 1024.0 * 3.0)))

        if fp >= 0.98:
            ret.append( (fp, patha, pathb) )

    if not ret:
        return None
    else:
        return ret


################################################################################
def compare_image_sims_pool(img_list):
    ret = []
    listing = []

    i = 1
    for simpath, path in img_list:
        listing.append( (simpath, path, img_list[i:]) )
        i += 1

    # Setup the pools?
    pool = Pool(processes=24)
    ret = pool.map(image_sim_compare_pool, listing)

    # Finish/close the pool
    pool.close()
    pool.join()

    # Flush out the None's
    ret = [item for item in ret if item != None]

    # Flatten the list
    ret = reduce(lambda x, y: x + y, ret)
    return ret


################################################################################
def cache_sim_data(sim):
    root = "/dev/shm/data"
    ret = []
    idx = 0

    for sim, pathb in sim:
        filename = os.path.join(root, str(idx) + ".sim")
        sim.dump(filename)
        idx += 1

        ret.append( (filename, pathb) )

    return ret


################################################################################
if __name__ == '__main__':
    # Option Parser
    usage = "usage: %prog [options] rootdir"
    parser = OptionParser(usage)
    parser.add_option("-e", "--exclude", dest="exclude", help="Exclude this directory")

    options, args = parser.parse_args()
    if len(args) != 1:
        parser.error("Need a root directory to recurse into")

    rootdir = args[0]
    exclude = options.exclude

    print "Generating image listing..."
    img = generate_img_list(rootdir, exclude)
    #calc_image_stats(img)

    # Generate SIM
#   h = hpy()
    print "Generating image sim..."
    sim = generate_sim_data(img)
#    print h.heap()

    print "Caching the sim data into a ramdisk..."
    sim_paths = cache_sim_data(sim)

    # Compare
#    h = hpy()
    print "Comparing image sim..."
    comp = compare_image_sims_pool(sim_paths)
#    print h.heap()

    # Flush out the 0.0's
    with open('out', 'w') as f:
        for val, path1, path2 in comp:
            print >>f, str(val) + " - " + path1 + " - " + path2
