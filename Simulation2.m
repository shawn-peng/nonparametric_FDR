function Simulation2(zeta1,zeta2)
%
n=10000;

x1 = zeta1.mu + zeta1.Delta*abs(randn(n,1)) + zeta1.Gamma*randn(n,1);
x2 = zeta2.mu + zeta2.Delta*abs(randn(n,1)) + zeta2.Gamma*randn(n,1);
subplot(3,1,1);
histogram(x1,'Normalization','pdf');
hold on;
histogram(x2,'Normalization','pdf');
hold off;
subplot(3,1,2);
xx = max(x1,x2);
r=rand(n,1);
ix = r >0.7;
yy= nan(n,1);
yy(ix)=x2(ix);
yy(~ix)=xx(~ix);

histogram(yy,'Normalization','pdf')
end