PRO bcd_phot_mosaic,mosaicfile,channel,outdir,sex=SEX,centroid=CENTROID,dumppos=DUMPPOS

; AUTHOR
;	John Livingston
; DATE
;	06/16/16

file_mkdir,outdir

;aper.pro setup
apr = [4, 6, 8, 10, 12]		;mosaic pixscale is 0.6", so this is a 2.4" aperture
skyrad = [24,40]
badpix = [0,0]

;conversion factors
ap_cor_ch1 = [1.208, 1.112, 1.070, 1.047, 1.032]
ap_cor_ch2 = [1.220, 1.112, 1.080, 1.048, 1.036]
;conv_fac = 35.174234		;MJy/sr --> uJy for native 1.2233"/pix
;conv_fac = 35.17515109		;MJy/sr --> uJy for 1.22304236283526 by 1.22361484682187
;conv_fac = 8.461933		;MJy/sr --> uJy for 0.6" mosaic pixels (taken from apex header)
conv_fac = 8.46159499		;MJy/sr --> uJy for 0.6" mosaic pixels (calculated)

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

if KEYWORD_SET(sex) then begin
	ss = strsplit(mosaicfile,'/',/extract)
	basename = strjoin(ss[0:n_elements(ss)-2],'/')
	radeclist = '/'+basename+'/radec.txt'
	;read in list of RA/Dec for sources to be measured in fitsfile image
	readcol,radeclist,ra,dec,format='D,D'
	;convert from sky to pixel coord.
	adxy,hdr,ra,dec,x,y
endif else begin
	hmin = 3
	fwhm = 3.33333
	roundlim = [-1.0,1.0]
	sharplim = [0.2,1.0]
	;find sources in mosaic image
	find, img, x, y, flux, sharp, roundness, hmin, fwhm, roundlim, sharplim
	;convert from pixel to sky coord.
	xyad,hdr,x,y,ra,dec
	;create ID
endelse

id = lindgen(n_elements(x))

;setup for writing results to disk

if KEYWORD_SET(dumppos) then begin
	pos_out = outdir+'/radec.txt'
	get_lun,pos
	openw,pos,pos_out
	for i=0,n_elements(x)-1 do begin

		if KEYWORD_SET(centroid) then begin
			; box_centroid,img,unc2,x[i],y[i],3,6,3,x0,y0,f0,b,xs,ys,fs,bs,c,cb,np
			gcntrd,img,x[i],y[i],x0,y0,3.33333,/silent
			; cntrd,img,x[i],y[i],x0,y0,3.33333,/silent
		endif else begin
			x0 = x[i]
			y0 = y[i]
		endelse

		xyad,hdr,x0,y0,ra0,dec0
		printf,pos,strtrim(strcompress(string([ra0,dec0])))

	endfor
	close,pos
endif else begin

	good_out = outdir+'/mosaic_phot_good.txt'
	bad_out = outdir+'/mosaic_phot_bad.txt'
	get_lun,good
	get_lun,bad
	openw,good,good_out,width=1200
	openw,bad,bad_out,width=1200

	printf,good,"id ra dec x y flux_mjy unc_mjy_aper unc_mjy_cbunc quality"
	printf,bad,"id ra dec x y"

	;loop through source pixel coordinates and do photometry at that location in image
	for i=0,n_elements(x)-1 do begin

		;check to make sure the pixel coordinates are finite, skip source if not
		if not finite(x[i]) or not finite(y[i]) then continue

		if KEYWORD_SET(centroid) then begin
			;centroid on x,y
			; box_centroid,img,unc2,x[i],y[i],3,6,3,x0,y0,f0,b,xs,ys,fs,bs,c,cb,np
			gcntrd,img,x[i],y[i],x0,y0,3.33333,/silent
			; cntrd,img,x[i],y[i],x0,y0,3.33333,/silent
		endif else begin
			x0 = x[i]
			y0 = y[i]
		endelse

		;get mean coverage per aperture
		aper,cov,x0,y0,cov_sum,cov_err,cov_sky,cov_skyerr,1,apr[0],badpix,$
			/flux,/nan,/exact,/silent,readnoise=0,setskyval=0,/meanback
		num_pix = 3.1415926 * apr[0] ^ 2
		mean_cov = cov_sum / num_pix

		;get photometry on centroid
		aper,img,x0,y0,flux_aper,fluxerr,sky,skyerr,phpadu*mean_cov,apr,$
			skyrad,badpix,/flux,/nan,/exact,/silent,readnoise=RONOISE/sqrt(mean_cov)

		;get uncertainties from unc mosaic
		aper,unc2,x0,y0,unc_sum,unc_err,unc_sky,unc_skyerr,1,apr[0],badpix,$
			/flux,/nan,/exact,/silent,readnoise=0,setskyval=0

		aper,unc2,x0,y0,unc_sum2,unc_err2,unc_sky2,unc_skyerr2,1,skyrad[0],badpix,$
			/flux,/nan,/exact,/silent,readnoise=0,setskyval=0

		aper,unc2,x0,y0,unc_sum3,unc_err3,unc_sky3,unc_skyerr3,1,skyrad[1],badpix,$
			/flux,/nan,/exact,/silent,readnoise=0,setskyval=0

		sky_area = 3.1415926 * (skyrad[1]^2 - skyrad[0]^2)
		sky_sum = unc_sum3 - unc_sum2
		sigma_tot = sqrt(unc_sum + sky_sum/sky_area)

		;convert flux and unc from MJy/sr to Jy and apply aperture correction
		flux_mjy = (flux_aper[0] * ap_cor[0] * conv_fac) * 1e-3
		unc_mjy = (fluxerr[0] * conv_fac) * 1e-3		; unc does not get aperture correction
		unc_mjy2 = (sigma_tot * conv_fac) * 1e-3

		;compute quality metric
		sum = 0
		for j=0,3 do sum += flux_aper[j] * ap_cor[j] / flux_aper[j+1] * ap_cor[j+1]
		quality = sum/4.

		;calculate RA/Dec of centroids
		xyad,hdr,x0,y0,ra0,dec0

		;print the data
		if finite(flux_mjy) eq 1 then begin
			printf,good,strtrim(strcompress([string(id[i]),$
				string([ra0,dec0,x0,y0,flux_mjy,unc_mjy,unc_mjy2,quality])]),1)
		endif else begin
			printf,bad,strtrim(strcompress([string(id[i]),$
				string([ra0,dec0,x0,y0])]),1)
		endelse

	endfor

	close,good
	close,bad

endelse

END
