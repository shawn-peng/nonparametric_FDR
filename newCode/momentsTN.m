function [m1,m2] = momentsTN(locs, scale)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here
x=locs/scale;
P=normcdf(x);
p=normpdf(x);
r=p./P;
r(P==0)=abs(x(P==0));
m1 = locs + scale*r;
m2 = locs.^2 + scale^2 + scale*locs.*r;
end

