PRO bcd_phot,cbcdfile,cbuncfile,maskfile,radeclist,channel,USE_MASK=use_mask

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
;	10/01/13
; MODIFIED
;	11/27/13 -- now supports new input format containing ID field per source

;bit flags to ignore
; ok_bitflags = [4, 16, 128, 20, 132, 144, 148]	;mask bits 2, 4, and 7
; ok_bitflags = [4, 16, 20]	;mask bits 2 and 4
ok_bitflags = [0]	;mask all bitflags

;aper.pro setup
apr = 3				;using BCDs so pixscale is native 1.2"/pix
; apr = 2				;using BCDs so pixscale is native 1.2"/pix
skyrad = [12,20]
badpix = [0,0]

;aperture photometry parameter string for pixel_phase_correct_gauss
ap_par = '3_12_20'
; ap_par = '2_12_20'

;conversion factors
ap_cor_ch1 = 1.112			;ch1 aperture correction 3 pix radius, 12-20 pix annulus
ap_cor_ch2 = 1.113			;ch2 aperture correction 3 pix radius, 12-20 pix annulus
; ap_cor_ch1 = 1.205			;ch1 aperture correction 2 pix radius, 12-20 pix annulus
; ap_cor_ch2 = 1.221			;ch2 aperture correction 2 pix radius, 12-20 pix annulus
conv_fac = 35.174234		;MJy/sr --> uJy for native 1.2233"/pix

if channel eq 1 then ap_cor = ap_cor_ch1
if channel eq 2 then ap_cor = ap_cor_ch2

;read in list of RA/Dec for sources to be measured in fitsfile image
readcol,radeclist,id,ra,dec,format='L,D,D'

;read in fits image and header
img = readfits(cbcdfile,hdr,/silent)

;read in uncertainty per pixel
unc = readfits(cbuncfile,/silent)

;read in imask file
if keyword_set(use_mask) then mask = readfits(maskfile,/silent)

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
if keyword_set(use_mask) then begin
	masked_out = work_dir+'masked_list.txt'
	get_lun,masked
	openw,masked,masked_out,width=1200,/append
endif

;loop through source pixel coordinates and do photometry at that location in image
for i=0,n_elements(x)-1 do begin

	;check to make sure the pixel coordinates are finite, skip source if not
	if not finite(x[i]) or not finite(y[i]) then continue

	;centroid on x,y
	box_centroider,img,unc^2,x[i],y[i],3,6,3,x0,y0,f0,b,xs,ys,fs,bs,c,cb,np

	;check for any flagged pixels inside the aperture, skip if so
	if keyword_set(use_mask) then begin
		ok_src = 0
		bitflag = badpix_aperture(mask,x0,y0,apr)
		if bitflag eq 0 then begin
			ok_src = 1
		endif else begin
			for j=0,n_elements(ok_bitflags)-1 do if bitflag eq ok_bitflags[j] then ok_src = 1
		endelse
		if not ok_src then begin
			printf,masked,strtrim(strcompress([string(id[i]),maskfile,$
				string([x0,y0]),string(bitflag)]),1)
			continue
		endif
	endif

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
	flux_uncor = flux_jy * 1e3

	;calculate RA/Dec of centroids
	xyad,hdr,x0,y0,ra0,dec0

	;print the data
	if finite(flux_aper) eq 1 then begin
		printf,good,strtrim(strcompress([string(id[i]),$
			string([ra[i],dec[i],ra0,dec0,x0,y0,flux_mjy,unc_mjy,flux_uncor])]),1)
	endif else begin
		printf,bad,strtrim(strcompress([string(id[i]),$
			string([ra[i],dec[i],ra0,dec0,x0,y0])]),1)
	endelse

endfor

close,good
close,bad
if keyword_set(use_mask) then close,masked

END