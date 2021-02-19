ds = {'S.cerevisiae_data.mat', 'S.cerevisiae3_data.mat', 'S.cerevisiae2_data.mat', ...
    'M.musculus_data.mat', 'M.musculus3_data.mat', 'M.musculus2_data.mat', 'H.sapiens3_data.mat', ...
    'H.sapiens2_data.mat','E.coli_data.mat', 'D.melanogaster_data.mat', 'A.thaliana_data.mat'};
ds = {'S.cerevisiae3_data.mat', 'S.cerevisiae2_data.mat', ...
    'M.musculus_data.mat', 'M.musculus3_data.mat', 'M.musculus2_data.mat', 'H.sapiens3_data.mat', ...
    'H.sapiens2_data.mat','E.coli_data.mat', 'D.melanogaster_data.mat', 'A.thaliana_data.mat'};
K=100;
for i=1:length(ds)
   S=load(['test_search/matdata/',ds{i}]);
    mat=S.mat2;
    if ~all(mat(:,1)==mat(:,2)) && length(mat)>100
        %[zeta,ll,lls,k] = EM12HNMax(mat,-1,1,K);
        S=load(['test_search/results/', ds{i}]);
        zeta=S.zeta; lls=S.lls; k=S.k;
        figure;
        subplot(1,2,1)
        plotFit(mat,zeta)
        subplot(1,2,2)
        [cdf,fdr,xm,t1p,fdr1p]=FDR(zeta);
        plot(xm,cdf,xm,fdr,'LineWidth',2)
        hold on
        xline(t1p,'-.',[num2str(fdr1p*100),'%','fdr:',num2str(t1p,4)])
        [cdfTrue,x]=sampleCDF(mat(:,1));
        plot(x,cdfTrue,'LineWidth',2,'LineStyle','--')
        cdfM=interp1(xm,cdf,x);
        cdfM=cdfM';
        cdfM=correctCDF(cdfM);
        deltaCDF=max(cdfM-cdfTrue);
        title(['\deltacdf:',num2str(deltaCDF,4)])
        hold off
        fn=ds{i};
        p=s1Dens(mat(:,1),zeta);
        ll1=mean(log(p));
        alpha=zeta.alpha;
        sgtitle({fn(1:length(fn)-4),['ll1: ',num2str(ll1)],['ll: ',num2str(lls(length(lls)))],['alpha:',num2str(alpha,3)]});
        saveas(gcf,['test_search/results/',fn(1:length(fn)-4),'.png'])
        %save(['test_search/results/', ds{i}],'zeta','llN','k','lls','ll1') 
    end
   
end