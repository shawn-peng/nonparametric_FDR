function [] = plotPdfs(s1,s2,p1,p1_c,p1_I1,p2,p2_c,p2_I1,p2_I2)
subplot(2,2,1); 
hold off
plot(s1,p1,s1,p1_c,s1,p1_I1,'LineWidth',1);
hold on;
histogram(s1,100,'Normalization','pdf');
hold off;
subplot(2,2,2);
hold off;
plot(s2,p2,s2,p2_c,s2,p2_I1,s2,p2_I2,'LineWidth',1);
hold on;
histogram(s2,100,'Normalization','pdf');
hold off;
end
