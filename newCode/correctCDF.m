function cdf = correctCDF(cdf)
cdf(cdf<0)=0;
cdf(cdf>1)=1;
il=find(cdf>0.9,1,'last');
n=length(cdf);
ix=1:n;
cdf(isnan(cdf)&ix>il)=1;
is=find(cdf<0.1,1,'first');
cdf(isnan(cdf)&ix<is)=0;
end

