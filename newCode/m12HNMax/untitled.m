ds = {'S.cerevisiae_data.mat', 'S.cerevisiae3_data.mat', 'S.cerevisiae2_data.mat', ...
    'M.musculus_data.mat', 'M.musculus3_data.mat', 'M.musculus2_data.mat', 'H.sapiens3_data.mat', ...
    'H.sapiens2_data.mat','E.coli_data.mat', 'D.melanogaster_data.mat', 'A.thaliana_data.mat'};
ds = {'S.cerevisiae3_data.mat', 'S.cerevisiae2_data.mat', ...
    'M.musculus_data.mat', 'M.musculus3_data.mat', 'M.musculus2_data.mat', 'H.sapiens3_data.mat', ...
    'H.sapiens2_data.mat','E.coli_data.mat', 'D.melanogaster_data.mat', 'A.thaliana_data.mat'};
ds = {'HeLa01ng_data.mat','HeLa1ng_data.mat', 'HeLa10ng_data.mat', 'HeLa50ng_data.mat','HeLa100ng_data.mat'};
%ds = {'c_elegans_data.mat','mouse_data.mat', 'human_data.mat', 'e_coli_data.mat','drosophila_data.mat'};
%ds = {'e_coli_data.mat'};
ds = {'HeLa01ng_data.mat','HeLa1ng_data.mat', 'HeLa10ng_data.mat', 'HeLa50ng_data.mat'};

B=50;


for dd=1:length(ds)
    d=ds{dd};
  S.ds=fn(1:(find(fn=='_')-1));
  S=load(['test_search/results/', ds{dd},'SNMax1','_B','.mat']);
  T(:,dd)=S.t1p;
end