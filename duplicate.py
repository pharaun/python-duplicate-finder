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
from optparse import OptionParser

################################################################################
import generate_sim
import compare_sim
import fhash
import generate_svd

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

    return total


################################################################################
def dup(comp, name):
    dup = calc_dup(comp)
    print "len " + name + ": " + str(len(comp)) + " dups: " + str(dup) + " nodup: " + str(len(comp) - dup)


################################################################################
def calc_dup(comp):
    # Assuming a common data storage format of this:
    # (file1, file2, (data... as needed))
    templist = list(comp)

    dupitem = len(comp) - len(set(comp))

    duppath = 0
    for patha, pathb, data in comp:
        if (patha == pathb):
            duppath += 1
            templist.remove( (patha, pathb, data) )

    dupinv = 0
    consumed_idx = []
    for idx, item in enumerate(templist):
        patha, pathb, data = item
        inv = (pathb, patha, data)

        if (inv in templist):
            if ((templist.index(inv)) not in consumed_idx):
                consumed_idx.append(idx)
                dupinv += 1

    return dupitem + duppath + dupinv


################################################################################
if __name__ == '__main__':
    # Option Parser
    usage = "usage: %prog [options] rootdir"
    parser = OptionParser(usage)
    parser.add_option("-e", "--exclude", dest="exclude", help="Exclude this directory")

    # Compare options
    parser.add_option("-s", "--sim", dest="imgsim", default=False, action='store_true',
            help="Compare with image similarity - abs(a-b)")
    parser.add_option("-f", "--filehash", dest="filehash", default=False, action='store_true',
            help="Compare with File hash")
    parser.add_option("-v", "--svd", dest="svd", default=False, action='store_true',
            help="Compare with Singular Value Decomposition")

    options, args = parser.parse_args()
    if len(args) != 1:
        parser.error("Need a root directory to recurse into")

    rootdir = args[0]
    imgsim = options.imgsim
    exclude = options.exclude
    filehash = options.filehash
    svd = options.svd

    if ((not imgsim) and (not filehash) and (not svd)):
        parser.error("Need to pick one compare option")

    print "Generating image listing..."
    start = time.time()
    img = generate_img_list(rootdir, exclude)
    print "Timing: " + str(time.time() - start) + " s"

    print
    print "Generating image stats..."
    total = calc_image_stats(img)

    comp4 = []
    if imgsim:
        # Generate SIM
        print
        print "Generating image sim..."
        start = time.time()
        sim = generate_sim.generate_sim_data(img)
        print "Timing: " + str(time.time() - start) + " s"

        # Diag output
        print
        print "Image stats..."
        print "Total images: " + str(total)
        print "Processed images: " + str(len(sim))
        print "Omitted images: " + str(total - len(sim))

        # Multiprocess Python Compare
        print
        print "Processing with c n-way compare..."
        start = time.time()
        comp4 = compare_sim.compare(sim)
        print "Timing: " + str(time.time() - start) + " s"

    elif filehash:
        # Generate file based hash & so forth for verification
        print
        print "Generating & Processing file based hashes..."
        start = time.time()
        comp4 = fhash.file_hashes(img)
        print "Timing: " + str(time.time() - start) + " s"

    elif svd:
        print
        print "Generating image svd..."
        start = time.time()
        sim = generate_svd.generate_svd_data(img)
        print "Timing: " + str(time.time() - start) + " s"

    # Detect duplicate.... duplicates
    print
    print "Detecting duplicate duplicates..."
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
