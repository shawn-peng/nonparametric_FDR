function d = areaBetweenCDF(cdf1,cdf2,s1)

dh = cdf1 - cdf2;
 
ddh = diff(dh);
ds = diff(s1);
n = size(dh,1);
s = dh(1:n-1) .* ds + 1/2 * ds .* ddh;
d = sum(abs(s));
end

