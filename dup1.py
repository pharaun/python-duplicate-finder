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
def generate_img_list(rootdir):
    img = []

    for root, sub, files in os.walk(rootdir):
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
            img_sum = []

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

                    # Append a rgb tuple to a list
                    img_sum.append( (r, g, b) )

                    i += x_inc
                j += y_inc

            return img_sum

    except IOError:
        return None


################################################################################
def generate_sim_data(img_list):
    img = []
#    total = len(img_list)
#    i = 0

    filepath = []
    for root, filename in img_list:
        filepath.append( os.path.join(root, filename))

    # Setup the pools?
    pool = Pool(processes=24)
    img = pool.map(generate_image_sim, filepath)

#    for root, filename in img_list:
#        img.append( (root, filename, generate_image_sim(os.path.join(root, filename))) )
#
#        i += 1
#        print "\r" + str(i) + "/" + str(total) + " - %.2f" % (100.0 * (float(i) / total)) + "%",
#        sys.stdout.flush()


    return img


if __name__ == '__main__':
    # Rest of the program
    if len(sys.argv) != 2:
        print sys.argv[0] + " <directory to recurse into>"
        sys.exit(2)
    rootdir = sys.argv[1]

    img = generate_img_list(rootdir)
    #calc_image_stats(img)

    h = hpy()
    sim = generate_sim_data(img)
    print h.heap()

    with open('out', 'w') as f:
        for line in sim:
            f.write(sim)


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
