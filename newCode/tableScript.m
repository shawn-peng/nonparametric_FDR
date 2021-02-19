ds = {'A.thaliana','D.melanogaster','E.coli','H.sapiens2', 'H.sapiens3', 'M.musculus',...
    'M.musculus2',  'M.musculus3','S.cerevisiae2', 'S.cerevisiae3'};
%ds = {'HeLa01ng_2','HeLa1ng'};

n=length(ds);
St=struct('GG',nan,'SNMix1',nan,'SNMix2',nan,'SNMax1',nan,'SNMax2',nan);
algo={'GG','SNMix1','SNMix2','SNMax1','SNMax2'};
algoAvail={'SNMax1','SNMax2'};
for a= 1:length(algo)
   t1p=[];t01p=[];t10p=[];deltaCdf=[];ll=[];
    for d =  1:length(ds)
        if any(strcmp(algo{a},algoAvail))
        S = jsondecode(fileread(['test_search/results/',ds{d},algo{a},'.json']));
        t1p=[t1p;S.t1p];
        deltaCdf=[deltaCdf;S.deltaCdf];
        ll=[ll;S.ll1];
        for ix=2:length(S.s1)
            S.fdr(ix)=max(min(S.fdr(ix),S.fdr(ix-1)),0);
        end
        [~,ii]=min(abs(S.fdr-0.001));
        t01p=[t01p;S.s1(ii)];
        [~,ii]=min(abs(S.fdr-0.1));
        t10p=[t10p;S.s1(ii)];
        else
            ll=ones(n,1)*-inf;
            deltaCdf=ones(n,1);
            t1p=ones(n,1)*inf;
            t01p=ones(n,1)*inf;
            t10p=ones(n,1)*inf;
        end
    end
    St.(algo{a})=struct('t01p',t01p,'t1p',t1p,'t10p',t10p,'deltaCdf',deltaCdf,'ll',ll);
end
save('test_search/results/tableFile.mat','St')