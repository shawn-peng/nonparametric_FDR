ds = {'S.cerevisiae_data.mat', 'S.cerevisiae3_data.mat', 'S.cerevisiae2_data.mat', ...
    'M.musculus_data.mat', 'M.musculus3_data.mat', 'M.musculus2_data.mat', 'H.sapiens3_data.mat', ...
    'H.sapiens2_data.mat','E.coli_data.mat', 'D.melanogaster_data.mat', 'A.thaliana_data.mat'};
ds = {'S.cerevisiae3_data.mat', 'S.cerevisiae2_data.mat', ...
    'M.musculus_data.mat', 'M.musculus3_data.mat', 'M.musculus2_data.mat', 'H.sapiens3_data.mat', ...
    'H.sapiens2_data.mat','E.coli_data.mat', 'D.melanogaster_data.mat', 'A.thaliana_data.mat'};
K=500;
B=3;
bp=0.1;
for i=1:length(ds)
   D=load(['test_search/matdata/',ds{i}]);
    mat=D.mat2;
    n=size(mat,1);
    S=struct();
    S.mat=mat;
    S.matb={};
    fn=ds{i}
    S.ds=fn(1:(find(fn=='_')-1));
    S.algo='SNMax2';
    if ~all(mat(:,1)==mat(:,2)) && length(mat)>100
        for i=1:B
            disp(i);
            matb=mat(randperm(n,ceil(bp*n)),:);
            S.matb{i}=matb;
            [S.zeta{i},S.ll(i),S.lls{i},S.k(i)] = EM23HNMax(matb,-1,1,K);
            %St=load(['test_search/results/', ds{i}]);
            %S=St.S;
            lls=S.lls{i};
            S.s1{i}=sort(matb(:,1));
            S.s2{i}=sort(matb(:,2));
            [S.pdfM{i},S.pdfC{i},S.pdfI1{i}]=s1Dens(S.s1{i},S.zeta{i});
            [S.fdr{i},S.cdf{i},xm,S.t1p(i),fdr1p(i)]=FDR(S.s1{i},S.zeta{i});
            [S.cdfTrue{i},x]=sampleCDF(S.s1{i});
            S.deltaCdf{i}=areaBetweenCDF(S.cdfTrue{i},S.cdf{i},S.s1{i});
            [S.pdfM_2,S.pdfC_2,S.pdfI1_1, S.pdfI1_2]=s2Dens(S.s2{i},S.zeta{i});
            
            S.ll1(i)=mean(log(S.pdfM{i}));
            fdr=S.fdr{i};
            S.pi_C(i)=1-fdr(1);
        
        end
            jsonStr = jsonencode(S);
            fid = fopen(['test_search/results/', S.ds ,S.algo,'_B','.json'], 'w');
            if fid == -1, error('Cannot create JSON file'); end
            fwrite(fid, jsonStr, 'char');
            fclose(fid);
            %saveas(gcf,['test_search/results/',fn(1:length(fn)-4),'.png'])
            save(['test_search/results/', S.ds,S.algo,'_B','.mat'],'S')
    end
end