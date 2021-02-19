St=load('test_search/results/tableFile.mat');
algo={'GG','SNMix1','SNMix2','SNMax1','SNMax2'};
%lAlgo={'GG','1SMix','2SMix','1SMax','2SMax'};
lAlgo={'\\GG','\\SNMixOne','\\SNMixTwo','\\SNMaxOne','\\SNMaxTwo'};
St=St.St;
algoAvail={'SNMax1','SNMax2'};
n=length(St.GG.ll);
t01p=nan(a,n);
t1p=nan(a,n);
t10p=nan(a,n);
deltaCdf=nan(a,n);
ll=nan(a,n);
aa=zeros(n,1);
for a= 1:length(algo)
   t01p(a,:)=St.(algo{a}).t01p';
   t1p(a,:)=St.(algo{a}).t1p';
   t10p(a,:)=St.(algo{a}).t10p';
   deltaCdf(a,:)=St.(algo{a}).deltaCdf';
   ll(a,:)=St.(algo{a}).ll';
end

[~,it01p]=min(t01p,[],1);
[~,it1p]=min(t1p,[],1);
[~,it10p]=min(t10p,[],1);
[~,idel]=min(deltaCdf,[],1);
[~,ill]=max(ll,[],1);
disp('\multirow{5}{5em}{0.1\% FDR Threshold}')
for a=1:length(algo)
    fprintf(['&',lAlgo{a}]);
    for i=1:n
        if it01p(i)==a
            fprintf(['& \\textbf{',num2str(t01p(a,i)),'}'])
        else
             fprintf(['&',num2str(t01p(a,i))])
        end
    end
   fprintf('\\\\ \n')
end
disp('\cline{1-12}')
disp('\multirow{5}{5em}{1\% FDR Threshold}')
for a=1:length(algo)
      fprintf(['&',lAlgo{a}]);
    for i=1:n
        if it1p(i)==a
            fprintf(['& \\textbf{',num2str(t1p(a,i)),'}'])
        else
             fprintf(['&',num2str(t1p(a,i))])
        end
    end
   fprintf('\\\\ \n')
end
disp('\cline{1-12}')
disp('\multirow{5}{5em}{10\% FDR Threshold}')
for a=1:length(algo)
    fprintf(['&',lAlgo{a}]);
    for i=1:n
        if it10p(i)==a
            fprintf(['& \\textbf{',num2str(t10p(a,i)),'}'])
        else
             fprintf(['&',num2str(t10p(a,i))])
        end
    end
   fprintf('\\\\ \n')
end
disp('\cline{1-12}')
disp('\multirow{5}{5em}{$\delta$ CDF}')
for a=1:length(algo)
     fprintf(['&',lAlgo{a}]);
    for i=1:n
        if idel(i)==a
            fprintf(['& \\textbf{',num2str(deltaCdf(a,i)),'}'])
        else
             fprintf(['&',num2str(deltaCdf(a,i))])
        end
    end
   fprintf('\\\\ \n')
end
disp('\cline{1-12}')
disp('\multirow{5}{5em}{log-likelihood}')
for a=1:length(algo)
   fprintf(['&',lAlgo{a}]);
    for i=1:n
        if ill(i)==a
            fprintf(['& \\textbf{',num2str(ll(a,i)),'}'])
        else
             fprintf(['&',num2str(ll(a,i))])
        end
    end
   fprintf('\\\\ \n')
end
disp('\cline{1-12}')