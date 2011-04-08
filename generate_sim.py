import os
import sys
import Image
import numpy as np
from multiprocessing import Pool

################################################################################
def generate_image_sim(path):
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
