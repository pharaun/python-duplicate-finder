import os
import sys
import hashlib

BLOCKSIZE = 1024 * 8
################################################################################
def file_hashes(img_list):
    img = []

    print "Chunking it into sizes..."
    # Break it up into "sizes" so we can save on disk reads latter for
    # the hashes
    sizedict = {}
    for root, filename in img_list:
        filepath = os.path.join(root, filename)
        size = os.path.getsize(filepath)

        files = sizedict.get(size, [])
        files.append(filepath)
        sizedict[size] = files

    # Cull the lists that are less than 2 item out of the sizedict
    for k, v in sizedict.items():
        if len(v) < 2:
            del sizedict[k]

    print "Hashing same size files..."
    # Perform the hashing on the remaining items
    hashdict = {}
    for file_list in sizedict.itervalues():
        for files in file_list:
            h = hashlib.md5() # Should look into other hashes, this will do and its quick

            with open(files, 'rb') as f:
                for chunk in iter(lambda: f.read(BLOCKSIZE), ''):
                    h.update(chunk)

            digest = h.hexdigest()
            hashes = hashdict.get(digest, [])
            hashes.append(files)
            hashdict[digest] = hashes

    # Cull the lists that are less than 2 item out of the hashdict
    for k, v in hashdict.items():
        if len(v) < 2:
            del hashdict[k]

    print "Comparing same hash files, byte by byte..."
    # Do byte to byte compare here to be 100% certain that two file that
    # hashes the same is really the same.
    to_process = []
    not_dupes = []

    for k, v in hashdict.items():
        for idxa in xrange(0, len(v)):
            for idxb in xrange ((idxa + 1), len(v)):
                tup = (v[idxa], v[idxb])
                tup2 = (v[idxb], v[idxa])

                if ((tup not in not_dupes) and (tup2 not in not_dupes)):
                    if ((tup not in img) and (tup2 not in img)):
                        to_process.append( tup )

        # Process the list of item to process
        for tup in to_process:
            fa, fb = tup

            try:
                f1 = open(fa, 'rb')
                f2 = open(fb, 'rb')

                while True:
                    buf1 = f1.read(BLOCKSIZE)
                    buf2 = f2.read(BLOCKSIZE)

                    if (buf1 and buf2):
                        # Compare
                        if (buf1 != buf2):
                            not_dupes.append( (fa, fb) )
                            break

                    elif (((not buf1) and buf2) or (buf1 and (not buf2))):
                        not_dupes.append( (fa, fb) )
                        break
                    else:
                        img.append( (fa, fb, None) )
                        break
            except:
                print "Unexpected error:", sys.exc_info()[0]
                raise
            finally:
                f1.close()
                f2.close()

        # Flush the to_process list
        to_process = []

    print "hash dup - " + str(len(hashdict))
    print "byte by byte dup - " + str(len(img))

    return img
