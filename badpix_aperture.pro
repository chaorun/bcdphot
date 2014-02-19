FUNCTION badpix_aperture,mask,x,y,radius

;checks the input mask (2d array) for the presence of
; nonzero pixels, returns 1 if any are within radius of
; (x,y)

xi = x-radius-1
xf = x+radius+1
yi = y-radius-1
yf = y+radius+1
dim = size(mask,/dim)
badpix = 0

for i=xi,xf do begin
	for j=yi,yf do begin
		if (i-x)^2 + (j-y)^2 lt radius^2 then begin
			if i ge 0 and j ge 0 and i le dim[0] and j le dim[1] then begin
				if mask[i,j] ne 0 then badpix = 1
			endif
		endif
	endfor
endfor

return,badpix

END