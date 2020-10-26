#!/usr/bin/python3
"""
	Code to implement Piantadosi 2018's "One parameter is always enough"
	This reads a png image with PIL and constructs a scatter plot by sampling y values at random from
	from the available x values. The code then encodes the y values into a single parameter theta, 
	following the bit-by-bit construction in the paper. 
"""

from random import sample, random
from PIL import Image

# Read an image
im = Image.open('miro.png')

Nx = 750 # How m any x-locations do we sample?

# construct the y points by sampling a black pixel for each 1:Nx spanning the range of the image
points = []
for x in range(Nx):

    px = int((float(x)/float(Nx)) * im.width) #NOTE: im.width may depend on the version of PIL?

    # find all of the y values that are black (or alomst black) at location px
    on = []
    for py in range(1,im.height-1):
        pixel = im.getpixel( (px,py) )
        if(pixel[-1] == 0): continue # transparent
        if all(z < 10 for z in pixel[:3]):
            on.append((im.height-float(py))/float(im.height))

        # fix missing with zero -- if there are no y values that are black, we pretend there was one at 0.01
        if(len(on) == 0):
            on.append( 0.01 )

        # choose which y value we sample
        yy = sample(on,1)[0]

        # append it to our list of points, scaling y (see paper) to be in (0,.5)
        points.append( (px, yy) )

# ------------------------------------------------------
# Define some useful functions, importing from mpmath
# ------------------------------------------------------

import mpmath
from mpmath import mp # mpmath allows arbitrary precision floating points

# set the precision required
mp.dps=2000

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
    return mpmath.asin(mpmath.sqrt(x)) / (2.0 * mpmath.pi)

def m(rx, theta):
    return mpmath.power(mpmath.sin(mpmath.power(2, rx) * mpmath.asin(mpmath.sqrt(theta))),2)

# ------------------------------------------------------
# Now the actual code
# ------------------------------------------------------

r = 8 # 2^{-r} is our precision

# concatenate the first r bits of each phiinv(y) to omaga represented as a string (omegastr)
omegastr = ""
for x,y in points:
    omegastr += float2binary(phiinv(y))[:r]

print( "# ", omegastr)
omega = binary2float(omegastr) # convert omega to a mp real 
theta = phi(omega) # compute theta
print( "# ", theta)

# now run the model to recover, print the fitted values
for x,y in points:
    ymodel = m(r*x,theta)

    # these fitted values may then be plotted
    print( x, y, float(ymodel) )
