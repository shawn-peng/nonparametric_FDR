function [] = plotPdfs(s1,p,p_c,p_I1)
figure;
subplot(1,2,1); 
hold off
plot(s1,p,s1,p_c,s1,p_I1,'LineWidth',2);
hold on;
histogram(s1,100,'Normalization','pdf');
hold off;
end
