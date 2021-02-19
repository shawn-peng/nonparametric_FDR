function mat = dataGen(alpha, beta, u_c, sigma_c, lambda_c, u_i1, sigma_i1, lambda_i1, u_i2, sigma_i2, lambda_i2, u_i3, sigma_i3, lambda_i3)
%
n=10000;
n1=int32(alpha*n);
n2=n-n1;
n3=int32(beta*n2);
n4=n2-n3;
delta_c= sqrt(lambda_c^2/(1+lambda_c^2));
delta_i1= sqrt(lambda_i1^2/(1+lambda_i1^2));
delta_i2= sqrt(lambda_i2^2/(1+lambda_i2^2));
delta_i3= sqrt(lambda_i3^2/(1+lambda_i3^2));

xc= u_c + sigma_c * delta_c * abs(randn(n,1)) + sigma_c*sqrt(1-delta_c^2)* randn(n,1);
xi1= u_i1 + sigma_i1 * delta_i1 * abs(randn(n,1)) + sigma_i1*sqrt(1-delta_i1^2)* randn(n,1);
xi2= u_i2 + sigma_i2 * delta_i2 * abs(randn(n,1)) + sigma_i2*sqrt(1-delta_i2^2)* randn(n,1);
xi3= u_i3 + sigma_i3 * delta_i3 * abs(randn(n,1)) + sigma_i3*sqrt(1-delta_i3^2)* randn(n,1);


s1 = [xc(1:n1);xi1(n1+1:n)];
s2 = [xi1(1:n1);xc(n1+1:n1+n3);xi2(n1+n3+1:n)];
s3 = [xi2(1:n1+n3);xi3(n1+n3+1:n)];
mat=[s1,s2,s3];
end

