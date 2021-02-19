ds = {'S.cerevisiae_data.mat', 'S.cerevisiae3_data.mat', 'S.cerevisiae2_data.mat', ...
    'M.musculus_data.mat', 'M.musculus3_data.mat', 'M.musculus2_data.mat', 'H.sapiens3_data.mat', ...
    'H.sapiens2_data.mat','E.coli_data.mat', 'D.melanogaster_data.mat', 'A.thaliana_data.mat'};
ds = {'S.cerevisiae3_data.mat', 'S.cerevisiae2_data.mat', ...
    'M.musculus_data.mat', 'M.musculus3_data.mat', 'M.musculus2_data.mat', 'H.sapiens3_data.mat', ...
    'H.sapiens2_data.mat','E.coli_data.mat', 'D.melanogaster_data.mat', 'A.thaliana_data.mat'};
ds = {'HeLa01ng_2_data.mat','HeLa1ng_data.mat'};
ds = {'M.musculus2_data.mat'};
ds = {'HeLa10ng_data.mat','HeL50ng_data.mat','HeL100ng_data.mat'};
K=500;
for i=1:length(ds)
   S=load(['test_search/matdata/',ds{i}]);
    mat=S.mat2;
    if ~all(mat(:,1)==mat(:,2)) && length(mat)>100
       [S.zeta,S.ll,S.lls,S.k] = EM23HNMax(mat,-1,-1,K);
       %S.zeta=zeta;S.ll=ll;S.lls=lls;S.k=k;
        fn=ds{i};
        fn=[fn(1:length(fn)-4),'2d'];
        %S=load(['test_search/results/', fn,'.mat']);
        %S=St.S;
        zeta=S.zeta; lls=S.lls; k=S.k;
        S.ll=S.lls(length(S.lls));
        S.mat=mat;
        S.s1=sort(mat(:,1));
        S.s2=sort(mat(:,2));
        [S.pdfM,S.pdfC,S.pdfI1]=s1Dens(S.s1,S.zeta);
        [S.fdr,S.cdf,xm,S.t1p,fdr1p]=FDR(S.s1,S.zeta);
        [S.cdfTrue,x]=sampleCDF(S.s1);
        S.deltaCdf=areaBetweenCDF(S.cdfTrue,S.cdf,S.s1);
        [S.pdfM_2,S.pdfC_2,S.pdfI1_2, S.pdfI2_2]=s2Dens(S.s2,S.zeta);
        S.ds=fn(1:(find(fn=='_')-1));
        S.algo='SNMax2';
        S.ll1=mean(log(S.pdfM));
        S.pi_C=1-S.fdr(1);
        jsonStr = jsonencode(S);
        fid = fopen(['test_search/results/', S.ds ,S.algo,'.json'], 'w');
        if fid == -1, error('Cannot create JSON file'); end
        fwrite(fid, jsonStr, 'char');
        fclose(fid);
        
        figure;
        subplot(2,2,1)
        plotFit(S.mat,S.zeta)
        subplot(2,2,3)
        plot(S.s1,S.cdf,S.s1,S.cdfTrue, S.s1,S.fdr,'LineWidth',2)
        hold on
        xline(S.t1p,'-.',[num2str(fdr1p*100),'%','fdr:',num2str(S.t1p,4)])
        title(['\deltacdf:',num2str(S.deltaCdf,4)])
        hold off
        saveas(gcf,['test_search/results/',fn,'.png'])
        save(['test_search/results/', fn,'.mat'],'S')
    end
   
end