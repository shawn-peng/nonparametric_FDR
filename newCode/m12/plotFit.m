function [] = plotFit(mat,zeta)
s1=sort(mat(:,1));
s1=s1(ix);
[p1,p1_c,p1_I1]=s1Dens(s1,zeta);
plot(s1,p1,s1,p1_c,s1,p1_I1,'LineWidth',2);
hold on;
histogram(s1,100,'Normalization','pdf');
hold off;
end

