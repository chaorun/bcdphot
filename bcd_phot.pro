PRO bcd_phot,cbcdfile,cbuncfile,radeclist,channel

; DESCRIPTION
;	Computes aperture photometry of the sources in radeclist
;	on the image in cbcdfile. Applies pixel phase correction but
;	not array location correction. Converts output to Jy and
;	prints a list of RA, Dec, x, y, flux, uncertainty to stdout.
; AUTHOR
;	John Livingston
; DATE
;	10/01/13

;aper.pro setup
;phpadu = 1          ;now using the correct value calculated from the header info
;apr = 2				;using BCDs so pixscale is native 1.2"/pix
apr = 3				;using BCDs so pixscale is native 1.2"/pix
skyrad = [12,20]
badpix = [0,0]

;aperture photometry parameter string for pixel_phase_correct_gauss
ap_par = '3_12_20'

;conversion factors
;ap_cor_ch1 = 1.205			;ch1 aperture correction 2 pix radius, 12-20 pix annulus
;ap_cor_ch2 = 1.221			;ch2 aperture correction 2 pix radius, 12-20 pix annulus
ap_cor_ch1 = 1.112			;ch1 aperture correction 3 pix radius, 12-20 pix annulus
ap_cor_ch2 = 1.113			;ch2 aperture correction 3 pix radius, 12-20 pix annulus
conv_fac = 33.847732		;MJy/sr --> uJy for native 1.2"/pix

if channel eq 1 then ap_cor = ap_cor_ch1
if channel eq 2 then ap_cor = ap_cor_ch2

;read in list of RA/Dec for sources to be measured in fitsfile image
readcol,radeclist,ra,dec,format='D,D'

;read in fits image and header
img = readfits(cbcdfile,hdr,/silent)

;read in uncertainty per pixel
unc = readfits(cbuncfile,/silent)

;calculate photons per digital unit from header values
GAIN = sxpar(hdr,'GAIN')
EXPTIME = sxpar(hdr,'EXPTIME')
FLUXCONV = sxpar(hdr,'FLUXCONV')
RONOISE = sxpar(hdr,'RONOISE')
phpadu = GAIN*EXPTIME/FLUXCONV

;convert from WCS to pixel coordinates
adxy,hdr,ra,dec,x,y

;loop through source pixel coordinates and compute photometry at that location in image
for i=0,n_elements(x)-1 do begin

	;centroid on x,y
	box_centroider,img,unc^2,x[i],y[i],3,6,3,x0,y0,f0,b,xs,ys,fs,bs,c,cb,np

	;get photometry on centroid
	aper,img,x0,y0,flux_aper,fluxerr,sky,skyerr,phpadu,apr,skyrad,badpix,/flux,/nan,/exact,/silent,readnoise=RONOISE

	;convert flux and unc from MJy/sr to Jy and apply aperture correction
	flux_jy = (flux_aper * ap_cor * conv_fac) * 1e-6
	unc = (fluxerr * ap_cor * conv_fac) * 1e-6

	;apply pixel phase correction
	corrected_flux = pixel_phase_correct_gauss(flux_jy,x0,y0,channel,ap_par)

	; ;fix annoying array dimension mismatch
	; flux = reform(flux)		
	; fluxerr = reform(fluxerr)	

	;calculate RA/Dec of centroids
	xyad,hdr,x0,y0,ra0,dec0

	;print the data
	print,strcompress([ra0,dec0,x0,y0,corrected_flux,unc])

endfor
	
END