
n = 10000;

alpha = 0.2;

zeta = {};

zeta.C.mu = 20;
zeta.C.sigma = 16;

zeta.I2.mu = -0.1;
zeta.I2.sigma = 5;

nz1 = int32(n*alpha);
nz0 = n - nz1;

xc = NSample(zeta.C, n);
xi2 = NSample(zeta.I2, n);
% xd = 

z = randi([0,1], n, 1);

thetaD2.mu = 0;
thetaD2.sigma = 5 + xi2/10;
xi1 = sampleHN(thetaD2, xi2, n);

s1 = xc .* z + xi1 .* (1-z);

% xi2 = xi1 - xd;
xd = xi1 - xi2;

s2 = xi1 .* z + xi2 .* (1-z);

s3 = xi2 .* z + (-1e8) .* (1-z);

mat = [s1, s2, s3];
[mat, I] = sort(mat, 2, 'descend');

s1 = mat(:,1);
s2 = mat(:,2);


figure
hold on;

histogram(s1, 50, 'Normalization', 'pdf')
histogram(s2, 50, 'Normalization', 'pdf')
legend(['s1'; 's2']);

figure
scatter(s1, s2)
xlabel('s1');
ylabel('s2');

figure
scatter(xi1, xd)
xlabel('i1');
ylabel('i1-i2');


figure
scatter(xi1, xi2)
xlabel('i1');
ylabel('i2');

figure
hold on;
histogram(xi1, 50, 'Normalization', 'pdf')
histogram(xi2, 50, 'Normalization', 'pdf')
legend(['i1'; 'i2']);

mat = [s1, s2];
