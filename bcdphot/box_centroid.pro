pro box_centroider, in_image, sigma2, xmax, ymax, boxwidth, backboxwidth, boxborder, $
                    x0, y0, f0, b, xs, ys, fs, bs, c, cb, np, MMM=mmm
;+
; NAME:
;    BOX_CENTROIDER
; 
; PURPOSE:
;    Calculates centroids for a source using a simple 1st moment box centroider
; 
; INPUTS:
;    IN_IMAGE: 2d float/double array containing input image, input image is presumed to 
;           be have bad pixels NaNed out
;    SIGMA2: 2d float/double array containing estimate of square of per pixel image 
;            uncertainty
;    XMAX: float/double scalar containing x position of peak pixel to centroid around
;    YMAX: float/double scalar containing y position of peak pixel to centroid around
;    BOXWIDTH: integer scalar, for size of centroiding box, actually halfboxwidth - 0.5
;    BACKBOXWIDTH: integer scalar, for size of rectangular background box
;    BOXBORDER: integer scalar, offset of background annulus in pixels from box aperture
;
; CONTROL KEYWORD:
;    MMM: if set, then backround estimate will be calculated using mmm.pro
;
; OUTPUTS:
;    X0: double scalar X position of centroid
;    Y0: double scalar Y position of centroid
;    F0: double scalar flux of source as summed in box after background removal
;    B: double scalar background level
;    XS: double scalar uncertainty in X centroid position
;    YS: double scalar uncertainty in Y centroid position
;    FS: double scalar uncertainty in measured flux accounting for uncertainty in background
;    BS: double scalar uncertainty in background
;    C: fixed scalar number of good pixels in aperture
;    CB: fixed scalar number of good pixels in background aperture
;    NP: fixed scalar number of noise pixels
;
; ALGORITHM:
;    
; HISTORY:
;    Recoding from get_centroids_nf procedure 26 May 2010 SJC
;    Revised input so that image is not directly operated on 24 Jan 2011 SJC
;    
;-

; Set sigma clipping threshold
	sigma = 3

; Copy input image to working array
	image = in_image

; Get image size
	sz = size(image)
	nx = sz[1]
	ny = sz[2]
		
; Arrays holding x and y coordinates of each pixel -- allow arbitrary sized 
; images
	xx = findgen(nx) # replicate(1.0, ny)
	yy = replicate(1.0, nx) # findgen(ny)
	
; Set box limits
	xa = xmax-boxwidth > 0
	xb = xmax+boxwidth < (nx-1)
	ya = ymax-boxwidth > 0
	yb = ymax+boxwidth < (ny-1)

; Set background regions
	xamin = xa - backboxwidth - boxborder > 0 < xa
	xamax = xa - boxborder > 0 < xa
	xbmax = xb + backboxwidth + boxborder > xb < (nx-1)
	xbmin = xb + boxborder > xb < (nx-1)
	yamin = ya - backboxwidth - boxborder > 0 < ya
	yamax = ya - boxborder > 0 < ya
	ybmax = yb + backboxwidth + boxborder > yb < (ny-1)
	ybmin = yb + boxborder > yb < (ny-1)

; Segment image into box and everything else
	mask = image * 0
	mask[xa:xb, ya:yb] = 1
	bptr = where(xx le xamax and xx ge xamin, bcount)
	if (bcount gt 0) then mask[bptr] = -1
	bptr = where(xx le xbmax and xx ge xbmin, bcount)
	if (bcount gt 0) then mask[bptr] = -1
	bptr = where(yy le yamax and yy ge yamin, bcount)
	if (bcount gt 0) then mask[bptr] = -1
	bptr = where(yy le ybmax and yy ge ybmin, bcount)
	if (bcount gt 0) then mask[bptr] = -1

; remove background - median of non box region
	bptr = where(mask eq -1 and image eq image, bcount)
	if (bcount gt 0) then begin 
		if (keyword_set(MMM)) then begin
			mmm, image[bptr], back, bsig
			bsig = bsig * bsig
		endif else begin
			back = median(image[bptr]) 
; One pass of outlier rejection for the background, could use mmm as well
			bmom = moment(image[bptr])
			bsig = bmom[1]
			bbptr = where(abs(image[bptr] - back) le sigma * sqrt(bsig), $
					               bbcount)
			if (bbcount gt 0) then begin
				back = median(image[bptr[bbptr]])
				bmom = moment(image[bptr[bbptr]] - back)
				bsig = bmom[1]
				bcount = bbcount
			endif else begin
				back = 0.0D
				bsig = 0.0D
			endelse
		endelse
	endif else begin
		back = 0.0D
		bsig = 0.0D
	endelse

; Perform background subtraction
	image = image - back
	b = back

; Perform centroid in 2*boxwidth+1 pixel box centered on peak pixel
	gptr = where(mask eq 1 and image eq image, gcount)
	fluxsum = total(image[gptr])
	flux2sum = total(image[gptr] * image[gptr])
	fluxsum2 = fluxsum * fluxsum
	x0 = total(xx[gptr] * float(image[gptr])) / fluxsum
	y0 = total(yy[gptr] * float(image[gptr])) / fluxsum
	
; Use summed flux as measure of total flux and calculate noise pixels		
	np = flux2sum / fluxsum
; Number of pixels in source and background regions
	c = float(gcount)
	cb = float(bcount)
	
; If a bad pixel exists in the source aperture return a flux of NaN.
	xptr = where(mask eq 1 and image ne image, xcount)
	if (xcount gt 0) then f = !VALUES.D_NAN else f = fluxsum

; Uncertainties
	xs = total((xx[gptr] - x0) * (xx[gptr] - x0) * sigma2[gptr])
	ys = total((yy[gptr] - y0) * (yy[gptr] - y0) * sigma2[gptr])
	xs = xs / fluxsum2
	ys = ys / fluxsum2
	bs = c * c * bsig / cb
	fs = total(sigma2[gptr]) + bs

return
end
