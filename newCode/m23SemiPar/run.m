ds = {'S.cerevisiae_data.mat', 'S.cerevisiae3_data.mat', 'S.cerevisiae2_data.mat', ...
    'M.musculus_data.mat', 'M.musculus3_data.mat', 'M.musculus2_data.mat', 'H.sapiens3_data.mat', ...
    'H.sapiens2_data.mat','E.coli_data.mat', 'D.melanogaster_data.mat', 'A.thaliana_data.mat'};
ds = {'S.cerevisiae3_data.mat', 'S.cerevisiae2_data.mat', ...
    'M.musculus_data.mat', 'M.musculus3_data.mat', 'M.musculus2_data.mat', 'H.sapiens3_data.mat', ...
    'H.sapiens2_data.mat','E.coli_data.mat', 'D.melanogaster_data.mat', 'A.thaliana_data.mat'};
res_dir = 'test_search/nonparam_results/';
addpath('newCode/')
addpath('newCode/m23SemiPar/')
addpath('newCode/distributions/')
addpath('newCode/NN/')
K=100;
% parpool(10);
parfor i=1:length(ds)
    S=load(['test_search/matdata_pride/',ds{i}]);
    mat=S.mat2;
    if ~all(mat(:,1)==mat(:,2)) && length(mat)>100
        figure;
%         [zeta,ll,lls,k] = EM23SPHNMax(mat,K,z);
        [zeta,ll,lls,k] = EM23SPHNMax_trueinit(mat,K,z,xc,xi1,xi2);
        fn=ds{i};
        fn=[fn(1:length(fn)-4),'2d'];
        %S=load(['test_search/results/', fn,'.mat']);
        %zeta=S.zeta; lls=S.lls; k=S.k;
        figure;
        binned_mat = zeta
        plotFit(mat,zeta)
        subplot(2,2,3)
        [fdr,cdf,xm,t1p,fdr1p]=FDR(mat(:,1),zeta);
        plot(xm,cdf,xm,fdr,'LineWidth',2)
        hold on
        xline(t1p,'-.',[num2str(fdr1p*100,1),'%','fdr:',num2str(t1p,'%.2f')])
        [cdfTrue,x]=sampleCDF(mat(:,1));
        plot(x,cdfTrue,'LineWidth',2,'LineStyle','--')
%         cdfM=interp1(xm,cdf,x);
%         cdfM=cdfM';
%         cdfM=correctCDF(cdfM);
%         cdfM=cdfM';
        cdfM=correctCDF(cdf);
        deltaCDF=max(cdfM-cdfTrue);
        title(['\delta cdf:',num2str(deltaCDF,4)])
        hold off
        p=s1Dens(mat(:,1),zeta);
        ll1=mean(log(p));
        alpha=zeta.alpha;
        sgtitle({fn(1:length(fn)-4),['ll1: ',num2str(ll1)],['ll: ',num2str(lls(length(lls)))],['alpha:',num2str(alpha,3)]});
        saveas(gcf,[res_dir,fn,'.png'])
%         save([res_dir, fn],'zeta','llN','k','lls','ll1') 
        parsave([res_dir, fn, '.mat'],zeta,ll,k,lls,ll1) 
    end
end
function parsave(fname,zeta,ll,k,lls,ll1)
  save(fname,'zeta','ll','k','lls','ll1')
end
