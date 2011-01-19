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
    sima, patha, idx = lista
    length = len(sim)
    ret = []

# This depends on "SIM" being available in the global space before
# forking
    while idx < length:
        simb, pathb = sim[idx]
#    for simb, pathb in listb:
        c = np.sum(np.absolute(np.subtract(sima, simb)))
        fp = (1.0 - (c / (255.0 * 1024.0 * 3.0)))

        if fp >= 0.98:
            ret.append( (fp, patha, pathb) )

# This depends on "SIM" being available in the global space before
# forking
        idx += 1

    if not ret:
        return None
    else:
        return ret


################################################################################
def compare_image_sims_pool(img_list):
    ret = []
    listing = []

    i = 1
    for arr1, path1 in img_list:
#        listing.append( (arr1, path1, img_list[i:]) )
# This depends on "SIM" being available in the global space before
# forking
        listing.append( (arr1, path1, i) )
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
    h = hpy()
    print "Generating image sim..."
    sim = generate_sim_data(img)
    print h.heap()


    # Compare
    h = hpy()
    print "Comparing image sim..."
    comp = compare_image_sims_pool(sim)
    print h.heap()

    # Flush out the 0.0's
    with open('out', 'w') as f:
        for val, path1, path2 in comp:
            print >>f, str(val) + " - " + path1 + " - " + path2



## Create a hash and compare the size
#size = {}
#
#for root, filename, what in img:
#    fs = os.path.getsize(os.path.join(root, filename))
#    if fs in size:
#        fslist = size[fs]
#        fslist.append((root, filename))
#        size[fs] = fslist
#    else:
#        size[fs] = [(root, filename)]
#
## Flush list in the dictionary that has only 1 file in it
#size = dict((k, v) for k, v in size.items() if len(v) > 1)
#
# Do hashes
#hashes = {}
#
#for root, filename, what in img:
#    sha512 = hashlib.sha512()
#
#    with open(os.path.join(root, filename)) as f:
#        for chunk in iter(lambda: f.read(128 * sha512.block_size), ''):
#            sha512.update(chunk)
#
#    digest = sha512.hexdigest()
#    if digest in hashes:
#        dlist = hashes[digest]
#        dlist.append((root, filename))
#        hashes[digest] = dlist
#    else:
#        hashes[digest] = [(root, filename)]
#
## Flush list in the dictionary that has only 1 file in it
#hashes = dict((k, v) for k, v in hashes.items() if len(v) > 1)
#
#################################################################################
#def image_sim_compare(sima, simb):
#    c = np.sum(np.absolute(np.subtract(sima, simb)))
#    return (1.0 - (c / (255.0 * 1024.0 * 3.0)))
#
#
#################################################################################
#def compare_image_sims(img_list):
#    ret = []
#
#    i = 1
#    for arr1, path1 in img_list:
#        for arr2, path2 in img_list[i:]:
#            ret.append( (image_sim_compare(arr1, arr2), path1, path2) )
#        i += 1
#
#    return ret
