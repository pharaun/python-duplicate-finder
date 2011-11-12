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
4) SVD
5) Various scaling algo, uses simple average value algo to scale down the
image for the SIM algo/feature

6) Image filtering ex:
    a) SIM does not deal good with images that is overall a solid block of
    single/majority color
    b) SIM also does not deal good with line drawing (because it'll average
	out to be white)
    c) Need to find a way to filter/identify these types of image to filter it
    out if possible before feeding it to various image processing algos


TODO:
1) Look in if i even need floats, may want to convert all data to uint8,
    including uint16 so then i can use the "MPSADBW" SSE4.1 instruction, which
    basically takes a whole bunch of integers and do "sum of absolute
    differences"

2) Look into uint8 loading instructions to get it loaded as fast as possible
_mm_sad_epu8


Clean up the #pragma, don't need the array or other stuff in there anymore so
makes sense to clean that stuff up


3) look into some sort of data structure that lets you identify/discard
already compared duplicate to cut down on compares a bit if possible
