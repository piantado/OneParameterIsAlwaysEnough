"""
	Code to implement Piantadosi 2018's "One parameter is always enough"
	This reads a png image with PIL and constructs a scatter plot by sampling y values at random from
	from the available x values. The code then encodes the y values into a single parameter theta, 
	following the bit-by-bit construction in the paper. 
    Modifications by Istvan Kolossvary May 2021, see comments.
"""
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("imageFile", help="image file")
parser.add_argument("onCutoff", help="cutoff to consider pixel lit, 0 (black) < cutoff < 1 (white), default = 0.1", nargs='?', default=0.1, type=float)
parser.add_argument("Nx", help="effective resolution in horizontal dimension, default = physical resolution", nargs='?', default=0, type=int)
parser.add_argument("r", help="2^{-r} is our precision in binary digits, default = 8", nargs='?', default=8, type=int)
parser.add_argument("sampleSize", help="size of random y samples in the vertical dimension for a given x, >= 0, default=0 means full sample", nargs='?', type=int)
args = parser.parse_args()

imageFile  = args.imageFile
onCutoff   = args.onCutoff
Nx         = args.Nx
r          = args.r
sampleSize = args.sampleSize

from random import sample, random
from PIL import Image

# Read an image
im = Image.open(imageFile)
width, height = im.size
print( "Image file           ", imageFile )
print( "Image width        = ", width )
print( "Image height       = ", height )
if( Nx == 0 ):
    Nx = width

# construct the y points by sampling a black pixel for each 1:Nx spanning the range of the image
points = []
count = 0
for x in range(Nx):

    px = int((float(x)/float(Nx)) * width) #NOTE: im.width may depend on the version of PIL?

    # find all of the y values that are black (or alomst black) at location px
    on = []
    for py in range(1,height-1):
        pixel = im.getpixel( (px,py) )
        if(pixel[-1] == 0): continue # transparent
        if all(z < onCutoff*255 for z in pixel[:3]):
            on.append((height-float(py))/float(height))

        """
        IK: I commented out the code block below, it has two major problems.
        1.  This block is incorrectly indented, it is part of the for loop above, which means that almost all 
            pixels will be processed, the vast majority of them empty, that takes an insane amount of time, 
            e.g., the Miro signature takes about 27 hours to compute on a fast CPU (and the result is wrong)
        2.  The original code can only handle a single y value per x pixel. This is why elephant.png and 
            miro.png draw their patterns with exteremely thick lines, so the algorithm can randomly pick a 
            single y value from multiple values that are "on" in a single x column (yy = sample(on,1)[0]). 
            Moreover, this algorithm cannot handle empty columns and must add a fake point at 0.01.
        
        The key to fix this is to introduce a counter as an additional element in the points tuple and for 
        computing 2^(rx) in the logistic map use this counter instead of px, which is still used for plotting and, then, 
        the theta function can handle multiple y values in a single column, plus there is no need to fudge empty columns.
        """

        # fix missing with zero -- if there are no y values that are black, we pretend there was one at 0.01
        #if(len(on) == 0):
            #on.append( 0.01 )
        # choose which y value we sample
        #yy = sample(on,1)[0]
        # append it to our list of points, scaling y (see paper) to be in (0,.5)
        #points.append( (px, yy) )

    if( sampleSize > 0 ):
        if(len(on) > 0 ):
            for _ in range( 1, min( sampleSize+1, len(on)+1)):
                yy = sample(on,1)[0]
                count += 1
                points.append( (count, px, yy) )
    else:
        for yOn in on:
            points.append( (count, px, yOn) )
            count += 1

print( "Nof points to fit  = ", len(points) )

# ------------------------------------------------------
# Define some useful functions, importing from mpmath
# ------------------------------------------------------

import mpmath
from mpmath import mp # mpmath allows arbitrary precision floating points
import math

# set the precision required
mp.prec = len(points) * r # binary digits

print( "Nof binary bits    = ", mp.prec )
print( "Nof decimal digits = ", mp.dps)

def float2binary(x):
    # convert x to a binary string
    assert 0 < x < 1

    outstr = ""
    while x>0:
        if x < 0.5:
            outstr = outstr+"0"
            x = 2*x
        else:
            outstr = outstr+"1"
            x = 2*x-1
    return outstr

def binary2float(s):
    # convert a binary string (after the decimal) to a float
    v = mp.mpf(0.0)
    for i, b in enumerate(s):
        if b == '1':
            v += mpmath.power(2,-(i+1))
    return v

def phi(z):
    return mpmath.power(mpmath.sin(z * mpmath.pi * 2.0), 2.0)

def phiinv(x):
    #return mpmath.asin(mpmath.sqrt(x)) / (2.0 * mpmath.pi)
    #return float(mpmath.asin(mpmath.sqrt(x)) / (2.0 * mpmath.pi))
    """
    IK: phiinv() is only used on float y values, there is no need for mpmath here!
        Moreover, float2binary() only needs float and not mpf, so returning mpf
        slows down the calculation in float2binary() prohibitively, AND the results are wrong.
    """
    return math.asin(math.sqrt(x)) / (2.0 * math.pi)

def m(rx, theta):
    #return mpmath.power(mpmath.sin(mpmath.power(2, rx) * mpmath.asin(mpmath.sqrt(theta))), 2)
    return mpmath.power(mpmath.sin(mpmath.power(2, rx) * m.asin_sqrt_theta), 2)


# ------------------------------------------------------
# Now the actual code
# ------------------------------------------------------

# make this a command line input parameter so it can be set 
# automatically based on the legth of the points array
#r = 8 # 2^{-r} is our precision

# concatenate the first r bits of each phiinv(y) to omega represented as a string (omegastr)
omegastr = ""
for n,x,y in points:
    omegastr += float2binary(phiinv(y))[:r]

print ( "Omega              = ", omegastr )
omega = binary2float(omegastr)  # convert omega to a mp real 
theta = phi(omega)              # compute theta
print ( "Theta              = ", theta )

m.asin_sqrt_theta = mpmath.asin(mpmath.sqrt(theta))  # IK: needs to be computed only once

# now run the model to recover, print the fitted values
for n,x,y in points:
    
    ymodel = m(r*n, theta)

    # these fitted values may then be plotted
    print (x, y, float(ymodel))
