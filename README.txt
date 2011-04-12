Basic requirements:
1) CPU that supports at least SSE3/SSSE3 for the hadd instruction for the SSE
enabled version, if not can remove the #define SSE out of c_sim.c file
2) Cython
3) Python 2.6.6 onward
4) Numpy
5) PIL

# Supported features:
1) File based hashing (md5) + byte by byte compare to verify the hash match
2) SIM (resample image to 32x32 then compute the difference) - if difference
is less than a certain amount they are classified as similiar

# Todo
1) SIFT/SUFT ?
2) various type of perceptual hashes
3) Wavlet....
