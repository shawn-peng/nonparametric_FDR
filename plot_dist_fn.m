function plot_dist_fn(mat, method, alpha, beta, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i, u_i2, sigma_i2, lambda_i2, u_i3, sigma_i3, lambda_i3)
minS = min(min(mat(mat~=-inf)));
maxS = max(max(mat));
x_values = minS:0.1:maxS;
omat = mat;
mat3 = mat;
mat2 = mat;
figure;
hold on;
y = skew_norm_pdf(x_values, u_c, sigma_c, lambda_c);
y2 = skew_norm_pdf(x_values, u_i, sigma_i, lambda_i);
y3 = alpha*y + (1-alpha)*y2;
% y4 = skew_norm_pdf
plot(x_values,y*alpha,'LineWidth',2);
plot(x_values,y2*(1-alpha),'LineWidth',2);
plot(x_values,y3,'LineWidth',2);

if (method == '_3s4c')
    y4 = skew_norm_pdf(x_values, u_i2, sigma_i2, lambda_i2);
    y5 = alpha*y2 + (1-alpha)*beta*y + (1-alpha)*(1-beta)*y4;
%     plot(x_values,y5,'LineWidth',2);
    y6 = skew_norm_pdf(x_values, u_i3, sigma_i3, lambda_i3);
    y7 = alpha*y4 + (1-alpha)*y6;
%     plot(x_values,y7,'LineWidth',2);
end

if (method == '_3s4c')
    histogram(mat3(:,1),300,'Normalization','pdf');
%     histogram(mat3(:,2),100,'Normalization','pdf');
%     histogram(mat3(:,3),100,'Normalization','pdf');
else
    histogram(omat(:,1),100,'Normalization','pdf');

    % histogram(mat(:,1),100,'Normalization','pdf');

    histogram(mat2(:,2),100,'Normalization','pdf');
end

legend({'dist\_correct'; 'dist\_incorrect'; 'mixture'; 'mixture2'; 'mixture3'; 'hist\_first'; 'hist\_second'});

hold off;
saveas(gcf,['test_search/distplot/',species,method,'.png'])

figure;
hold on;

S1 = sort(omat(:,1));
m = size(S1,1);
h1 = (1:m) / m;
pemp = plot(S1,h1);

h11 = skew_norm_cdf(x_values, u_c, sigma_c, lambda_c);
h12 = skew_norm_cdf(x_values, u_i, sigma_i, lambda_i);
h1e = alpha*h11 + (1-alpha)*h12;
pest = plot(x_values, h1e);

seq = flip(unique(omat(:,1)));

n = size(seq,1);
for i = 1:n
    s = seq(i);
    fdr = fdr_x(alpha, beta, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i, s);
%     fdr
    if fdr > 0.1
        break
    end
end
xline(s);
text(s,0.05,'\leftarrow 1%fdr')
xlabel('-log(EValue)');
ylabel('CDF');
ylim([0,1]);

legend([pemp, pest], {'empirical CDF'; 'estimated CDF'});

saveas(gcf,['test_search/fitting curve/',species,method,'.png'])

end
