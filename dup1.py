#!/usr/bin/env python
"""
Python Duplicate picture finder, this mainly tests
various stuff related to the duplicate finding, such as
file/directory access, imhdr(identification) and so forth
"""

import os
import sys
import time
import Image
import imghdr
import hashlib
import numpy as np
from multiprocessing import Pool
from optparse import OptionParser

################################################################################
import compare_sim

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
def dup(comp, name):
    dup = calc_dup(comp)
    print "len " + name + ": " + str(len(comp)) + " dups: " + str(dup) + " nodup: " + str(len(comp) - dup)


################################################################################
def calc_dup(comp):
    templist = list(comp)

    dupitem = len(comp) - len(set(comp))

    duppath = 0
    for fp, patha, pathb in comp:
        if (patha == pathb):
            duppath += 1
            templist.remove( (fp, patha, pathb) )

    dupinv = 0
    consumed_idx = []
    for idx, item in enumerate(templist):
        fp, patha, pathb = item
        inv = (fp, pathb, patha)

        if (inv in templist):
            if ((templist.index(inv)) in consumed_idx):
                #skip
                1+1
            else:
                consumed_idx.append(idx)
                dupinv += 1

    return dupitem + duppath + dupinv


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
    start = time.time()
    img = generate_img_list(rootdir, exclude)
    print "Timing: " + str(time.time() - start) + " s"
    print "Generating image stats..."
    print
    calc_image_stats(img)
    print

    # Generate SIM
    print "Generating image sim..."
    start = time.time()
    sim = generate_sim_data(img)
    print "Timing: " + str(time.time() - start) + " s"

    # Multiprocess Python Compare
    print "Processing with c n-way compare..."
    start = time.time()
    comp4 = compare_sim.compare(sim)
    print "Timing: " + str(time.time() - start) + " s"

    # Detect duplicate.... duplicates
    print
    dup(comp4, "c")

#    for idx in xrange(0,len(comp1)):
#        fpa, pathaa, pathba = comp1[idx]
#        fpb, pathab, pathbb = comp2[idx]
#
#        if (fpa != fpb) or (pathaa != pathab) or (pathba != pathbb):
#            print fpa, fpb, " - ", pathaa, pathab, " - ", pathba, pathbb
#
#    # Flush out the list
#    with open('comp1', 'w') as f:
#        for val, path1, path2 in comp1:
#            print >>f, str(val) + " - " + path1 + " - " + path2
#
#    # Flush out the list
#    with open('comp2', 'w') as f:
#        for val, path1, path2 in comp2:
#            print >>f, str(val) + " - " + path1 + " - " + path2
