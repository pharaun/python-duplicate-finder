import os
import sys
import Image, ImageOps
import numpy as np
from multiprocessing import Pool

################################################################################
def generate_image_simaaaa(path):
    mode = ('RGB',)
    conv = ('1', 'P', 'L', 'LA', 'RGBA', 'RGBX', 'CMYK',)

    try:
        with open(path) as f:
            image = Image.open(f)
            if not (image.mode in mode):
                if not (image.mode in conv):
                    return None
                else:
                    # First check to see if its an animated image if so
                    # omit it - This only works for animated GIF's...
                    if image.info.has_key('version'):
                        if image.info['version'].__contains__('GIF'):
                            if image.info.has_key('duration'):
                                if image.info['duration'] > 0:
                                    return None

                    # Convert it - Ideally better to process each on
                    # their own but this will work for extending
                    # coverage to other image types
                    try:
                        image = image.convert('RGB')
                    except:
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
    except:
        print "Unexpected error:", sys.exc_info()[0]
        raise


################################################################################
def generate_sim_dataaaaa(img_list):
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
def generate_svd_data(img_list):
    mode = ('RGB',)

    filepath = []
    for root, filename in img_list:
        filepath.append( os.path.join(root, filename))

    try:
        for path in filepath:
            with open(path) as f:
                image = Image.open(f)
                if not (image.mode in mode):
                    return None

                grey = ImageOps.grayscale(image)
                grey2 = np.asarray(grey.resize((210, 315), Image.NEAREST))
                s1 = np.linalg.svd(grey2, compute_uv=False)

                grey2 = np.asarray(grey.resize((210, 315), Image.BILINEAR))
                s2 = np.linalg.svd(grey2, compute_uv=False)

                grey2 = np.asarray(grey.resize((210, 315), Image.BICUBIC))
                s3 = np.linalg.svd(grey2, compute_uv=False)

                grey2 = np.asarray(grey.resize((210, 315), Image.ANTIALIAS))
                s4 = np.linalg.svd(grey2, compute_uv=False)

                print np.sum(s1)
                print np.sum(s2)
                print np.sum(s3)
                print np.sum(s4)

                print np.sum(np.fabs(s1 - s2))
                print np.sum(np.fabs(s1 - s3))
                print np.sum(np.fabs(s1 - s4))
                print np.sum(np.fabs(s2 - s3))
                print np.sum(np.fabs(s2 - s4))
                print np.sum(np.fabs(s3 - s4))

                print (np.sum(np.fabs(s1 - s2)) / (210 * 315 * 255))*100.0
                print (np.sum(np.fabs(s1 - s3)) / (210 * 315 * 255))*100.0
                print (np.sum(np.fabs(s1 - s4)) / (210 * 315 * 255))*100.0
                print (np.sum(np.fabs(s2 - s3)) / (210 * 315 * 255))*100.0
                print (np.sum(np.fabs(s2 - s4)) / (210 * 315 * 255))*100.0
                print (np.sum(np.fabs(s3 - s4)) / (210 * 315 * 255))*100.0

                print np.allclose(s1, s2, atol=150)
                print np.allclose(s1, s3, atol=150)
                print np.allclose(s1, s4, atol=100)
                print np.allclose(s2, s3, atol=100)
                print np.allclose(s2, s4, atol=100)
                print np.allclose(s3, s4, atol=100)

                return None

                # Do stuff
                (np.rint(np.nan_to_num(image_pixel)), path)

    except IOError:
        return None
    except:
        print "Unexpected error:", sys.exc_info()[0]
        raise
