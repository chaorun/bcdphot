PRO bcd_phot_mosaic,mosaicfile,channel,sextractor=SEXTRACTOR

; DESCRIPTION
;	Computes aperture photometry of the sources in radeclist
;	on the image in cbcdfile. Applies pixel phase correction but
;	not array location correction. Converts output to Jy and
;	writes RA, Dec, x, y, flux, uncertainty, and uncorrected flux
;	to disk, as well as details of any sources for which aper.pro
;	returned a NaN -- either due to the source being off the array
;	or the presence of bad pixels.
; AUTHOR
;	John Livingston
; DATE
;	05/26/16

;bit flags to ignore
; ok_bitflags = [4, 16, 128, 20, 132, 144, 148]	;mask bits 2, 4, and 7
; ok_bitflags = [4, 16, 20]	;mask bits 2 and 4
; ok_bitflags = [0]	;mask all bitflags

;aper.pro setup
; apr = 3				;using BCDs so pixscale is native 1.2"/pix
; apr = 2				;using BCDs so pixscale is native 1.2"/pix
apr = 4					;mosaic pixscale is 0.6", so this is a 2.4" aperture
; skyrad = [12,20]
skyrad = [24,40]
badpix = [0,0]

;aperture photometry parameter string for pixel_phase_correct_gauss
; ap_par = '3_12_20'
; ap_par = '2_12_20'

;conversion factors
; ap_cor_ch1 = 1.112			;ch1 aperture correction 3 pix radius, 12-20 pix annulus
; ap_cor_ch2 = 1.113			;ch2 aperture correction 3 pix radius, 12-20 pix annulus
ap_cor_ch1 = 1.205			;ch1 aperture correction 2 pix radius, 12-20 pix annulus
ap_cor_ch2 = 1.221			;ch2 aperture correction 2 pix radius, 12-20 pix annulus
;conv_fac = 35.174234		;MJy/sr --> uJy for native 1.2233"/pix
conv_fac = 8.461933			;MJy/sr --> uJy for mosaic pixels

if channel eq 1 then ap_cor = ap_cor_ch1
if channel eq 2 then ap_cor = ap_cor_ch2

;read in fits image and header
img = readfits(mosaicfile,hdr,/silent)

;read in uncertainty per pixel
;unc = readfits(uncfile,/silent)

;calculate photons per digital unit from header values
GAIN = sxpar(hdr,'GAIN')
EXPTIME = sxpar(hdr,'EXPTIME')
FLUXCONV = sxpar(hdr,'FLUXCONV')
RONOISE = sxpar(hdr,'RONOISE')
phpadu = GAIN*EXPTIME/FLUXCONV


if KEYWORD_SET(sextractor) then begin
	ss = strsplit(mosaicfile,'/',/extract)
	basename = strjoin(ss[0:n_elements(ss)-2],'/')
	radeclist = '/'+basename+'/radec.txt'
	;read in list of RA/Dec for sources to be measured in fitsfile image
	readcol,radeclist,ra,dec,format='D,D'
	;convert from sky to pixel coord.
	adxy,hdr,ra,dec,x,y
endif else begin
	hmin = 3
	fwhm = 4
	roundlim = [-1.0,1.0]
	sharplim = [0.2,1.0]
	;find sources in mosaic image
	find, img, x, y, flux, sharp, roundness, hmin, fwhm, roundlim, sharplim
	;convert from pixel to sky coord.
	xyad,hdr,x,y,ra,dec
	;create ID
endelse

id = indgen(n_elements(x))

;setup for writing results to disk
; sep = strsplit(radeclist,'/')
; work_dir = strmid(radeclist,0,sep[n_elements(sep)-1])
work_dir = './'
good_out = work_dir+'good_list.txt'
bad_out = work_dir+'bad_list.txt'
masked_out = work_dir+'masked_list.txt'
get_lun,good
get_lun,bad
get_lun,masked
openw,good,good_out,width=1200
openw,bad,bad_out,width=1200
openw,masked,masked_out,width=1200

;loop through source pixel coordinates and do photometry at that location in image
for i=0,n_elements(x)-1 do begin

	;check to make sure the pixel coordinates are finite, skip source if not
	if not finite(x[i]) or not finite(y[i]) then continue

	x0 = x[i]
	y0 = y[i]

	;get photometry on centroid
	aper,img,x0,y0,flux_aper,fluxerr,sky,skyerr,phpadu,apr,skyrad,badpix,$
		/flux,/nan,/exact,/silent,readnoise=RONOISE

	;convert flux and unc from MJy/sr to Jy and apply aperture correction
	flux_jy = (flux_aper * ap_cor * conv_fac) * 1e-6
	unc_jy = (fluxerr * conv_fac) * 1e-6		; unc does not get aperture correction

	;apply pixel phase correction
	flux_mjy = flux_jy * 1e3
	unc_mjy = unc_jy * 1e3

	;calculate RA/Dec of centroids
	xyad,hdr,x0,y0,ra0,dec0

	;print the data
	if finite(flux_aper) eq 1 then begin
		printf,good,strtrim(strcompress([string(id[i]),$
			string([ra[i],dec[i],ra0,dec0,x0,y0,flux_mjy,unc_mjy])]),1)
	endif else begin
		printf,bad,strtrim(strcompress([string(id[i]),$
			string([ra[i],dec[i],ra0,dec0,x0,y0])]),1)
	endelse

endfor

close,good
close,bad
close,masked

END
