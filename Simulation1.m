function Simulation1(zeta1,zeta2)
%
n=10000;

x1 = zeta1.mu + zeta1.Delta*abs(randn(n,1)) + zeta1.Gamma*randn(n,1);
x2 = zeta2.mu + zeta2.Delta*abs(randn(n,1)) + zeta2.Gamma*randn(n,1);
subplot(3,1,1)
histogram(x1,'Normalization','pdf')
hold on;
histogram(x2,'Normalization','pdf')
hold off;
subplot(3,1,2)
xx = x1-x2;
xxp = xx(xx>0);
xxn = abs(xx(xx<0));
histogram(xxp,'Normalization','pdf')
hold on;
histogram(xxn,'Normalization','pdf')
hold off;
end