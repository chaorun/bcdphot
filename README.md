bcdphot
=======

High precision BCD photometry pipeline for large Spitzer surveys, 
optimized for use on multicore machines. 

Applies the corrections for the pixel phase effect as well as the
array location-dependent photometric effect, both of which must be
applied at the individual BCD level instead of on a single mosaic.
See [here](http://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/warmimgcharacteristics/) for details.

Supports IRAC's HDR mode, and utilizes the Spitzer Science Center
pipeline imask files to discard contaminated measurements.

Requires [IDL](http://www.exelisvis.com/ProductsServices/IDL.aspx)
and [astrolib](http://idlastro.gsfc.nasa.gov/).

Python requirements:
- numpy
- scipy
- statsmodels
- multiprocessing
- simplejson
- pyyaml
- pyfits
- pywcs
