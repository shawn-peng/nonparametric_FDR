ds = {'S.cerevisiae_data.mat', 'S.cerevisiae3_data.mat', 'S.cerevisiae2_data.mat', ...
    'M.musculus_data.mat', 'M.musculus3_data.mat', 'M.musculus2_data.mat', 'H.sapiens3_data.mat', ...
    'H.sapiens2_data.mat','E.coli_data.mat', 'D.melanogaster_data.mat', 'A.thaliana_data.mat'};
 K=1000;
for i=1:length(ds)
   S=load(['test_search/matdata/',ds{i}]);
    mat=S.mat2;
    if ~all(mat(:,1)==mat(:,2)) && length(mat)>100
        [zeta,ll,lls,k] = EM23HNMax(mat,-1,1,K);
        figure;
        plotFit(mat,zeta)
        fn=ds{i};
        p=s1Dens(mat(:,1),zeta);
        ll1=mean(log(p));
        sgtitle({fn(1:length(fn)-4),['ll1:',num2str(ll1)],['ll:',num2str(ll)]});
        saveas(gcf,['test_search/results/',fn(1:length(fn)-4),'.png'])
        save(['test_search/results/', ds{i}],'zeta','llN','k','lls','ll1') 
    end
   
end