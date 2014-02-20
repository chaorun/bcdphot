FUNCTION get_bits,in

in = long(in)
out = bytarr(32,n_elements(in))
for i=0,31 do out(31-i,*)=(in and 2L^i)/2L^i
bits = 31-where(out eq 1)
return,reverse(bits)

END