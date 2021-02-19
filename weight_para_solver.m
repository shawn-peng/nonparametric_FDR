syms alpha beta gamma delta epsilon eta1 eta2 eta3
syms C1 I11 C2 I21 I22 I32 I33
syms M1 M2 M3

% Ab = [
%     %alpha, beta,   gamma,  delta,  epsi,   et1,    eta2,   eta3,   b
%     C1+I21, 0,      0,      0,      0,      1,      0,      -1,     0;
%     0,      C2,     0,      0,      0,      1,      -1,     -1,     0;
%     0,      0,      I22+I33,0,      0,      1,      -1,     0,      0;
%     0,      0,      0,      I11,    0,      0,      1,      0,      0;
%     0,      0,      0,      0,      I32,    0,      0,      1,      0;
%     1,      1,      1,      0,      0,      0,      0,      0,      1;
%     1,      1,      1,      0,      0,      0,      0,      0,      1;
%     1,      1,      1,      0,      0,      0,      0,      0,      1;
%     ]

I11 = M1 - C1;
I22 = M2 - C2 - I21;
I33 = M3 - I32;
% assume(C1 + I11 == M1);
% assume(C2 + I21 + I22 == M2);
% assume(I32 + I33 == M3);


eq1 = 1/alpha * (C1 + I21) + eta1 - eta3 == 0;
eq2 = 1/beta * C2 + eta1 - eta2 - eta3 == 0;
eq3 = 1/gamma * (I22 + I33) + eta1 - eta2 == 0;
eq4 = 1/delta * I11 + eta2 == 0;
eq5 = 1/epsilon * I32 + eta3 == 0;
eq6 = alpha + beta + gamma == 1;
eq7 = beta + gamma == delta;
eq8 = alpha + beta == epsilon;

Y = solve(eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8, ...
    alpha, beta, gamma, delta, epsilon, eta1, eta2, eta3);
