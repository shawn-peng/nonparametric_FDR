function [p1,p1_c,p1_I1,p2,p2_c,p2_I1,p2_I2] = s1s2Dens(s1,s2,zeta)

n = 10000;
xc = sample(zeta.C, n);
xi2 = sample(zeta.I2, n);

z = rand(n, 1) < zeta.alpha;
theta.mu = zeta.D1.mu;
theta.sigma = sqrt(zeta.D1.ss);
xi1 = sampleHN(theta, xi2, n);

xs1 = xc .* z + xi1 .* (1-z);

% xi2 = xi1 - xd;
xd = xi1 - xi2;

xs2 = xi1 .* z + xi2 .* (1-z);

xs3 = xi2 .* z + (-1e8) .* (1-z);

mat = [xs1, xs2, xs3];
[mat, I] = sort(mat, 2, 'descend');

xs1 = mat(:,1)';
xs2 = mat(:,2)';


z = logical(z);
xc1 = xc(z & I(:,1)==1);
x1i1 = xi1(~z | (z & I(:,1)~=1));

xc2 = xc(z & I(:,1)==2);
x2i1 = xi1(z & I(:,1)==1);
x2i2 = xi2((z & I(:,3)==1) | (~z & I(:,1)==1));

alpha = sum(z & I(:,1)==1) / n;
beta = sum(z & I(:,1)==2) / n;

figure;
hold off;
h = histogram(xs1,50,'Normalization','pdf');
hold on;
[freq_count, bin_edges] = histcounts(x1i1,h.BinEdges,'Normalization','pdf');
bin_values = (bin_edges(1:end-1) + bin_edges(2:end))/2;
% b = bar(bin_values,freq_count*(1-alpha),1);
% b.FaceAlpha = 0.5;
plot(bin_values,freq_count*(1-alpha),'LineWidth',2);
[freq_count, bin_edges] = histcounts(xc1,h.BinEdges,'Normalization','pdf');
bin_values = (bin_edges(1:end-1) + bin_edges(2:end))/2;
% b = bar(bin_values,freq_count*alpha,1);
% b.FaceAlpha = 0.5;
plot(bin_values,freq_count*(alpha),'LineWidth',2);
% histogram(xc1,50,'Normalization','pdf','BinEdges',h.BinEdges);
legend({'s1'; 'i1'; 'c'});



figure;
hold off;
h = histogram(xs2,50,'Normalization','pdf');
hold on;
[freq_count, bin_edges] = histcounts(x2i1,h.BinEdges,'Normalization','pdf');
bin_values = (bin_edges(1:end-1) + bin_edges(2:end))/2;
% b = bar(bin_values,freq_count*(1-alpha),1);
% b.FaceAlpha = 0.5;
plot(bin_values,freq_count*(alpha),'LineWidth',2);
[freq_count, bin_edges] = histcounts(x2i2,h.BinEdges,'Normalization','pdf');
bin_values = (bin_edges(1:end-1) + bin_edges(2:end))/2;
% b = bar(bin_values,freq_count*(1-alpha),1);
% b.FaceAlpha = 0.5;
plot(bin_values,freq_count*(1-alpha-beta),'LineWidth',2);
[freq_count, bin_edges] = histcounts(xc2,h.BinEdges,'Normalization','pdf');
bin_values = (bin_edges(1:end-1) + bin_edges(2:end))/2;
% b = bar(bin_values,freq_count*alpha,1);
% b.FaceAlpha = 0.5;
plot(bin_values,freq_count*(beta),'LineWidth',2);
% histogram(xc1,50,'Normalization','pdf','BinEdges',h.BinEdges);
legend({'s2'; 'i1'; 'i2'; 'c'});



figure;
hold off;
h = histogram(s1,50,'Normalization','pdf');
hold on;
[freq_count, bin_edges] = histcounts(x1i1,h.BinEdges,'Normalization','pdf');
bin_values = (bin_edges(1:end-1) + bin_edges(2:end))/2;
% b = bar(bin_values,freq_count*(1-alpha),1);
% b.FaceAlpha = 0.5;
plot(bin_values,freq_count*(1-alpha),'LineWidth',2);
[freq_count, bin_edges] = histcounts(xc1,h.BinEdges,'Normalization','pdf');
bin_values = (bin_edges(1:end-1) + bin_edges(2:end))/2;
% b = bar(bin_values,freq_count*alpha,1);
% b.FaceAlpha = 0.5;
plot(bin_values,freq_count*(alpha),'LineWidth',2);
% histogram(xc1,50,'Normalization','pdf','BinEdges',h.BinEdges);
legend({'s1'; 'i1'; 'c'});



figure;
hold off;
h = histogram(s2,50,'Normalization','pdf');
hold on;
[freq_count, bin_edges] = histcounts(x2i1,h.BinEdges,'Normalization','pdf');
bin_values = (bin_edges(1:end-1) + bin_edges(2:end))/2;
% b = bar(bin_values,freq_count*(1-alpha),1);
% b.FaceAlpha = 0.5;
plot(bin_values,freq_count*(alpha),'LineWidth',2);
[freq_count, bin_edges] = histcounts(x2i2,h.BinEdges,'Normalization','pdf');
bin_values = (bin_edges(1:end-1) + bin_edges(2:end))/2;
% b = bar(bin_values,freq_count*(1-alpha),1);
% b.FaceAlpha = 0.5;
plot(bin_values,freq_count*(1-alpha-beta),'LineWidth',2);
[freq_count, bin_edges] = histcounts(xc2,h.BinEdges,'Normalization','pdf');
bin_values = (bin_edges(1:end-1) + bin_edges(2:end))/2;
% b = bar(bin_values,freq_count*alpha,1);
% b.FaceAlpha = 0.5;
plot(bin_values,freq_count*(beta),'LineWidth',2);
% histogram(xc1,50,'Normalization','pdf','BinEdges',h.BinEdges);
legend({'s2'; 'i1'; 'i2'; 'c'});






[freq_count, bin_value] = hist(xc1,50);


end

