function [m1, m2] = trunc_norm_moments(U, sigma)
%trunc_norm_pdf pdf of truncated norm distribution
%   Detailed explanation goes here
%     cdf = normcdf(U/sigma);
%     flags = cdf == 0;
%     cdf(flags) = min(min(cdf(~flags)));
%     pdf = normpdf(U/sigma);
%     m1 = U + sigma * (normpdf(U/sigma) ./ cdf);
%     m2 = U.^2 + sigma^2 + sigma * U .* (normpdf(U/sigma) ./ cdf);
%     m2(isnan(m2)) = inf;
[m1,m2]=momentsTN(U,sigma);
end

