FUNCTION get_bits,in

; some of this code was scavenged from:
; http://www.astro.washington.edu/docs/idl/cgi-bin/getpro/library32.html?DEC2BIN

in = long(in)
out = bytarr(32,n_elements(in))
for i=0,31 do out(31-i,*)=(in and 2L^i)/2L^i
bits = 31-where(out eq 1)
return,reverse(bits)

END