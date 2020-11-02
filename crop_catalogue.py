#!/usr/bin/env python

from astropy import wcs
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np

from beam_value_at_radec import beam_value
from beam_value_at_radec import parse_metafits

from astropy.io import fits
from optparse import OptionParser

usage="Usage: %prog [options] <file>\n"
parser = OptionParser(usage=usage)
parser.add_option('--ra',type="string", dest="ra",
                    help="The RA center of your crop in decimal degrees")
parser.add_option('--dec',type="string", dest="dec",
                    help="The declination center of your crop in decimal degrees")
parser.add_option('--radius',type="float", dest="radius",
                    help="The radius of your crop circle in degrees")
parser.add_option('--minflux',type="float", dest="minflux", default=1.0,
                    help="The minimum flux density of the source (default = 1.0 Jy)")
parser.add_option('--attenuate', action="store_true", dest="attenuate", default=False,
                    help="Use the primary beam to attenuate the catalogue fluxes when calculating the minimum flux (default=False); NB: this will be slow for a very low flux density cut-off!")
parser.add_option('--metafits',type="string", dest="metafits",
                    help="The metafits file for your observation (only necessary if using the beam attenuation option)")
parser.add_option('--catalogue',type="string", dest="cat",
                    help="The filename of the catalogue you want to read in (fits format).")
parser.add_option('--racol',type="string", dest="racol", default="RAJ2000",
                    help="The name of the RA column (in decimal degrees) in your catalogue (default = RAJ2000)")
parser.add_option('--decol',type="string", dest="decol", default="DEJ2000",
                    help="The name of your Dec column (in decimal degrees) in your catalogue (default = DEJ2000)")
parser.add_option('--fluxcol',type="string", dest="fluxcol", default="int_flux_wide",
                    help="The name of the flux density column (default = int_flux_wide)")
parser.add_option('--alphacol',type="string", dest="alphacol", default="alpha",
                    help="The name of the spectral index column (default = alpha)")
parser.add_option('--output',type="string", dest="output", default="cropped_catalogue.fits",
                    help="The filename of the output cropped catalogue (fits format).")
(options, args) = parser.parse_args()

pos = SkyCoord(options.ra, options.dec, unit=(u.deg, u.deg))

# Read in the catalogue
temp = fits.open(options.cat)
data = temp[1].data

coords = SkyCoord(data[options.racol], data[options.decol], unit=(u.deg, u.deg))
separations = np.array(pos.separation(coords))
# Select sources within crop radius of your original RA and Dec
indices = np.where(separations<options.radius)
# Select only sources which meet the minimum flux density criterion
fluxcut = np.where(data[options.fluxcol]>options.minflux)
# Select only sources with non-zero spectral indices
alphacut = np.where(np.logical_not(np.isnan(data[options.alphacol])))
indices = np.intersect1d(fluxcut,indices)
indices = np.intersect1d(alphacut,indices)

# Now that we have some subset, check if beam attenuation is on, and if it is, attenuate and recalculate:

if options.attenuate:
    if options.metafits:
        t, delays, freq = parse_metafits(options.metafits)
        x,y = beam_value(data[options.racol][indices], data[options.decol][indices], t, delays, freq)
# This might be more correct for Stokes values, but I'm not sure it's what I want...
#        i = (x**2 + y**2)/2
# Also, at some point, I would like to save the full Jones matrix into the source model
        i = (x + y)/2
# Attenuate the fluxes
        data[options.fluxcol][indices] = i * data[options.fluxcol][indices]
# Select only the sources which meet the minimum flux criterion
        subindices = np.where(data[options.fluxcol][indices] > options.minflux)
        indices = indices[subindices]
    else:
        print("No metafits file selected.")
        sys.exit(1)

nselected=indices.shape[0]
noriginal=data.shape[0]

if indices.shape > 0:
    # Write out the sources
    temp[1].data = data[indices]
    temp.writeto(options.output,overwrite=True)
    print("Selected {0} of {1} sources".format(nselected,noriginal))
else:
    print "No sources selected!"
