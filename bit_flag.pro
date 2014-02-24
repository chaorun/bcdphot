FUNCTION bit_flag,bits

sum = 0
for i=0,n_elements(bits)-1 do sum += 2^bits[i]
return, sum

END