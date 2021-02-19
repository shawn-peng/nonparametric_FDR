
clear

% list_species = {
%     'H.sapiens3';
%     'M.musculus2';
%     'M.musculus3';
%     };
% clear

list_species = {
%     'C.elegans2';
    'S.cerevisiae';
    'M.musculus2';
    'M.musculus3';
    };

% species = 'M.musculus'
% species = 'H.sapiens'
% species = 'H.sapiens2'
% species = 'H.sapiens3'
% species = 'H.sapiens4'
% species = 'C.elegans'
% species = 'D.melanogaster'
% species = 'S.cerevisiae'
% species = 'S.cerevisiae2'
% species = 'S.cerevisiae3'
% species = 'E.coli'
% species = 'A.thaliana'

% species

n = size(list_species, 1)
for i = 1:n
    species = list_species{i};
    run parse_psm.m
    [alpha, beta, u_c, sigma_c, lambda_c, u_i, sigma_i, lambda_i] = EM2_2(omat',1,1)
    run plot_dist.m
    run plot_fdr.m
end
