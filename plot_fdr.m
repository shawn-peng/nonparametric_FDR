

curv = [0,0,0];

seq = flip(unique(omat(:,1)));

n = size(seq,1);
for i = 1:n
    s = seq(i);
    fdr = fdr_x(alpha, beta, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i, s);
%     fdr
    if fdr > 0.1
        break
    end
    curv(i,:) = [fdr, sum(omat(:,1)>s), s];
end

figure();
plot(curv(:,1), curv(:,2));

% csvwrite('twomix.csv', curv);
dlmwrite(['test_search/',species,method,'.csv'], curv, 'delimiter', ',', 'precision', '%25.20f');
