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
; MODIFIED
;	11/27/13 -- now supports new input format containing ID field per source

;aper.pro setup
apr = 3				;using BCDs so pixscale is native 1.2"/pix
skyrad = [12,20]
badpix = [0,0]

;aperture photometry parameter string for pixel_phase_correct_gauss
ap_par = '3_12_20'

;conversion factors
ap_cor_ch1 = 1.112			;ch1 aperture correction 3 pix radius, 12-20 pix annulus
ap_cor_ch2 = 1.113			;ch2 aperture correction 3 pix radius, 12-20 pix annulus
; conv_fac = 33.847732		;MJy/sr --> uJy for native 1.2"/pix
; conv_fac = 35.17421909111017		;MJy/sr --> uJy for native 1.2233"/pix
conv_fac = 35.174234		;MJy/sr --> uJy for native 1.2233"/pix

if channel eq 1 then ap_cor = ap_cor_ch1
if channel eq 2 then ap_cor = ap_cor_ch2

;read in list of RA/Dec for sources to be measured in fitsfile image
readcol,radeclist,id,ra,dec,format='L,D,D'

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

;setup for writing results to disk
sep = strsplit(radeclist,'/')
work_dir = strmid(radeclist,0,sep[n_elements(sep)-1])
good_out = work_dir+'good_list.txt'
bad_out = work_dir+'bad_list.txt'
get_lun,good
get_lun,bad
openw,good,good_out,width=1200,/append
openw,bad,bad_out,width=1200,/append

;loop through source pixel coordinates and do photometry at that location in image
for i=0,n_elements(x)-1 do begin

	;check to make sure the pixel coordinates are finite, skip source if not
	if not finite(x[0]) or not finite(y[0]) then continue

	;centroid on x,y
	box_centroider,img,unc^2,x[i],y[i],3,6,3,x0,y0,f0,b,xs,ys,fs,bs,c,cb,np

	;get photometry on centroid
	aper,img,x0,y0,flux_aper,fluxerr,sky,skyerr,phpadu,apr,skyrad,badpix,$
		/flux,/nan,/exact,/silent,readnoise=RONOISE

	;convert flux and unc from MJy/sr to Jy and apply aperture correction
	flux_jy = (flux_aper * ap_cor * conv_fac) * 1e-6
	unc = (fluxerr * ap_cor * conv_fac) * 1e-6

	;apply pixel phase correction
	corrected_flux = pixel_phase_correct_gauss(flux_jy,x0,y0,channel,ap_par)
	flux_mjy = corrected_flux * 1e3
	unc_mjy = unc * 1e3

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

END