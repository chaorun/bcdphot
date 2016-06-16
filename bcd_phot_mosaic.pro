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
;conv_fac = 35.17515109		;MJy/sr --> uJy for 1.22304236283526 by 1.22361484682187
;conv_fac = 8.461933			;MJy/sr --> uJy for 0.6" mosaic pixels (taken from apex header)
conv_fac = 8.46159499			;MJy/sr --> uJy for 0.6" mosaic pixels (calculated)

if channel eq 1 then ap_cor = ap_cor_ch1
if channel eq 2 then ap_cor = ap_cor_ch2

;read in fits image and header
img = readfits(mosaicfile,hdr,/silent)

;read in uncertainty per pixel
ss = strsplit(mosaicfile,'/',/extract)
basename = strjoin(ss[0:n_elements(ss)-2],'/')
uncfile = '/'+basename+'/mosaic_unc.fits'
unc = readfits(uncfile,/silent)
unc2 = unc^2 		;square the uncertainties for quadrature sum

;read in the coverage map
ss = strsplit(mosaicfile,'/',/extract)
basename = strjoin(ss[0:n_elements(ss)-2],'/')
covfile = '/'+basename+'/mosaic_cov.fits'
cov = readfits(covfile,/silent)

;calculate photons per digital unit from header values
; GAIN = sxpar(hdr,'GAIN')
; EXPTIME = sxpar(hdr,'EXPTIME')
; FLUXCONV = sxpar(hdr,'FLUXCONV')
; RONOISE = sxpar(hdr,'RONOISE')
; phpadu = GAIN*EXPTIME/FLUXCONV
RONOISE = 15.0		;per BCD
phpadu = 306.126	;per BCD
;use the coverage map to compute the average coverage within each
;aperture, then adjust the above values appropriately, i.e. the
;read noise decreases by sqrt(n) but the exposure time increases by n


if KEYWORD_SET(sextractor) then begin
	ss = strsplit(mosaicfile,'/',/extract)
	basename = strjoin(ss[0:n_elements(ss)-2],'/')
	radeclist = '/'+basename+'/radec.txt'
	;read in list of RA/Dec for sources to be measured in fitsfile image
	readcol,radeclist,ra,dec,format='D,D'
	;convert from sky to pixel coord.
	adxy,hdr,ra,dec,x,y
endif else begin
	hmin = 5
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
get_lun,good
get_lun,bad
openw,good,good_out,width=1200
openw,bad,bad_out,width=1200

;loop through source pixel coordinates and do photometry at that location in image
for i=0,n_elements(x)-1 do begin

	;check to make sure the pixel coordinates are finite, skip source if not
	if not finite(x[i]) or not finite(y[i]) then continue

	x0 = x[i]
	y0 = y[i]

	;get mean coverage per aperture
	aper,cov,x0,y0,cov_sum,cov_err,cov_sky,cov_skyerr,1,apr,badpix,$
		/flux,/nan,/exact,/silent,readnoise=0,setskyval=0
	num_pix = 3.1415926 * apr ^ 2
	mean_cov = cov_sum / num_pix

	;get photometry on centroid
	aper,img,x0,y0,flux_aper,fluxerr,sky,skyerr,phpadu*mean_cov,apr,$
		skyrad,badpix,/flux,/nan,/exact,/silent,readnoise=RONOISE/sqrt(mean_cov)

	;get uncertainties from unc mosaic
	aper,unc2,x0,y0,unc_sum,unc_err,unc_sky,unc_skyerr,1,apr,badpix,$
		/flux,/nan,/exact,/silent,readnoise=0,setskyval=0

	;convert flux and unc from MJy/sr to Jy and apply aperture correction
	flux_mjy = (flux_aper * ap_cor * conv_fac) * 1e-3
	unc_mjy = (fluxerr * conv_fac) * 1e-3		; unc does not get aperture correction
	unc_mjy2 = (sqrt(unc_sum) * conv_fac) * 1e-3

	; flux_mjy = flux_jy * 1e3
	; unc_mjy = unc_jy * 1e3
	; unc_mjy2 = unc_jy2 * 1e3

	;calculate RA/Dec of centroids
	xyad,hdr,x0,y0,ra0,dec0

	;print the data
	if finite(flux_aper) eq 1 then begin
		printf,good,strtrim(strcompress([string(id[i]),$
			string([ra[i],dec[i],ra0,dec0,x0,y0,flux_mjy,unc_mjy,unc_mjy2])]),1)
	endif else begin
		printf,bad,strtrim(strcompress([string(id[i]),$
			string([ra[i],dec[i],ra0,dec0,x0,y0])]),1)
	endelse

endfor

close,good
close,bad

END
