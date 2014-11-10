PRO xsc_phot,fitsfile,radeclist,exp

img = readfits(fitsfile, hdr)

readcol,radeclist,id,ra,dec,r_ext,format='A,D,D,D'

if exp eq 'long' then begin
	phpadu = 696.8874700718277
	ronoise = 7.8
endif 
if exp eq 'short' then begin
	phpadu = 23.62330407023145
	ronoise = 22.6
endif
badpix = [0,0]

adxy,hdr,ra,dec,x,y
radii = r_ext / sxpar(hdr,'PXSCAL2')

conv_fac = 8.461933

for i=0,n_elements(x)-1 do begin

	skyrad = [radii[i] * 1.5, radii[i] * 2.0]
;	print,'using:',radii[i],skyrad[0],skyrad[1]
	aper,img,x[i],y[i],flux_aper,fluxerr,sky,skyerr,phpadu,radii[i],$
		skyrad,badpix,/flux,/nan,/exact,/silent,readnoise=RONOISE
	flux_mjy = flux_aper * conv_fac * 1e-3
	unc_mjy = fluxerr * conv_fac * 1e-3
	sky_mjy = sky * conv_fac * 1e-3
	skyerr_mjy = skyerr * conv_fac * 1e-3

	print,strtrim( strcompress( $
		; [id[i], string([ra[i],dec[i],x[i],y[i],flux_mjy,unc_mjy,sky_mjy,skyerr_mjy])] ) )
		[id[i], string([flux_mjy,unc_mjy,sky_mjy,skyerr_mjy])] ) )

endfor

END