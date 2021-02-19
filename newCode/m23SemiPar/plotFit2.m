function [] = plotFit(mat,zeta)
s1=sort(mat(:,1));
s2=sort(mat(:,2));
% s3=sort(mat(:,3));
% [p1,p1_c,p1_I1]=s1Dens(s1,zeta);
% [p2,p2_c,p2_I1,p2_I2]=s2Dens(s2,zeta);
% %[p3,p3_I2,p3_I3]=s3Dens(s3,zeta);
% subplot(2,2,1); 
% hold off
% plot(s1,p1,s1,p1_c,s1,p1_I1,'LineWidth',2);
% hold on;
% histogram(s1,100,'Normalization','pdf');
% hold off;
% subplot(2,2,2);
% hold off;
% plot(s2,p2,s2,p2_c,s2,p2_I1,s2,p2_I2,'LineWidth',2);
% hold on;
% histogram(s2,100,'Normalization','pdf');
% hold off;

% subplot(2,2,3); 
% plot(s3,p3,s3,p3_I2,s3,p3_I3,'LineWidth',2);
% hold on;
% histogram(s3,100,'Normalization','pdf');
% hold off;

% [p1,p1_c,p1_I1]=s1Dens(s1,zeta);
% [p2,p2_c,p2_I1,p2_I2]=s2Dens(s2,zeta);

[p1,p1_c,p1_I1,p2,p2_c,p2_I1,p2_I2] = s1s2Dens(s1, s2, zeta);

subplot(2,2,3);
hold off
plot(s1,p1,s1,p1_c,s1,p1_I1,'LineWidth',2);
hold on;
histogram(s1,100,'Normalization','pdf');
hold off;
legend('mix', 'C', 'I1', 's1');
subplot(2,2,4);
hold off;
plot(s2,p2,s2,p2_c,s2,p2_I1,s2,p2_I2,'LineWidth',2);
hold on;
histogram(s2,100,'Normalization','pdf');
hold off;
legend('mix', 'C', 'I1', 'I2', 's2');
end

