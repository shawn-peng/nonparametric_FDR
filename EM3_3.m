function [alpha, beta, u_c, sigma_c, lambda_c, u_i1, sigma_i1, lambda_i1, u_i2, sigma_i2, lambda_i2, u_i3, sigma_i3, lambda_i3] = EM3_3(S,slc,sl1,sl2,sl3)
%EM The EM algorithm to estimate parameters
%   consider the first and the second score are dependent, they can not be
%   true at the same time
% tollerance = 1e-8;
tollerance = 1e-7;
[N, M] = size(S);

ll = 0;
prev_ll = -1;

S1_sorted = sort(S(1,:), 'descend');
S2_sorted = sort(S(2,:), 'descend');
S3_sorted = sort(S(3,:), 'descend');
% S_sorted = S;

prior_thres = 30;
nc1 = sum(S1_sorted > prior_thres);
nc2 = sum(S2_sorted > prior_thres);
% alpha = 0.1;
% alpha = nc1 / M;
alpha = 0.5;
% beta = 0.05;
beta = (nc2 / M) / (1 - alpha);
beta=0.5;

[u_c, sigma_c, lambda_c] = sn_para_est(S1_sorted(1:int32(M*alpha)));
[u_i1, sigma_i1, lambda_i1] = sn_para_est(S1_sorted(int32(M*alpha):M));
[u_i2, sigma_i2, lambda_i2] = sn_para_est(reshape(S2_sorted, 1, []));
[u_i3, sigma_i3, lambda_i3] = sn_para_est(reshape(S3_sorted, 1, []));
lambda_c = lambda_c * slc;
lambda_i1 = lambda_i1 * sl1;
lambda_i2 = lambda_i2 * sl2;
lambda_i3 = lambda_i3 * sl3;
%lambda_c=0;
%lambda_i1=0;
%lambda_i2=0;
%lambda_i3=0;
%sigma_c=1;
%sigma_i1=1;
%sigma_i2=1;
%sigma_i3=1;
%u_c=10;
%u_i1=0;
%u_i2=-10;
%u_i3=-20;
% sigma_c = sigma_i2;
% u_i3 = u_i2;
% sigma_i3 = sigma_i2;
% lambda_i3 = abs(lambda_i2) * sl3;

u_i = u_i1;
sigma_i = sigma_i1;
lambda_i = lambda_i1;

s1 = S(1,:);
%s1 = s1(s1~=0);
s2 = S(2,:);
%s2 = s2(s2~=0);
s3 = S(3,:);
%s3 = s3(s3~=0);

M1 = size(s1,2);
M2 = size(s2,2);
M3 = size(s3,2);
% u_c = 5;
% sigma_c = 1;
% lambda_c = 4;
% u_i = 0;
% sigma_i = 1;
% lambda_i = 2;
% alpha = 0.5;

while abs(ll - prev_ll) > tollerance
    prev_ll = ll;
    
%     P_c = skew_norm_pdf(S, u_c, sigma_c, lambda_c);
%     P_i = skew_norm_pdf(S, u_i, sigma_i, lambda_i);
    
    % inner sum sum across rows (k = 1, N) producing 1 x M vector
%     ll = sum( sum(log(P_i), 1) + log(alpha * sum(P_c ./ P_i, 1) + (1-alpha) * N) );
    ll = func_ll3_3(s1, s2, s3, alpha, beta, u_c, sigma_c, lambda_c, u_i1, sigma_i1, lambda_i1, u_i2, sigma_i2, lambda_i2, u_i3, sigma_i3, lambda_i3);
    %disp(ll);
    disp(ll - prev_ll);
    %disp([alpha, beta, u_c, sigma_c, lambda_c, u_i1, sigma_i1, lambda_i1, u_i2, sigma_i2, lambda_i2, u_i3, sigma_i3, lambda_i3]);
    
%     plot_dist_fn(S, '_3s4c', alpha, beta, u_c, sigma_c, lambda_c, u_i1, sigma_i1, lambda_i1, u_i2, sigma_i2, lambda_i2, u_i3, sigma_i3, lambda_i3);
    
    delta_i1 = lambda_i1 / sqrt(1+lambda_i1^2);
    Delta_i1 = sigma_i1 * delta_i1;
    Gamma_i1 = sigma_i1^2 * (1-delta_i1^2);
    delta_i2 = lambda_i2 / sqrt(1+lambda_i2^2);
    Delta_i2 = sigma_i2 * delta_i2;
    Gamma_i2 = sigma_i2^2 * (1-delta_i2^2);
    delta_i3 = lambda_i3 / sqrt(1+lambda_i3^2);
    Delta_i3 = sigma_i3 * delta_i3;
    Gamma_i3 = sigma_i3^2 * (1-delta_i3^2);
    delta_c = lambda_c / sqrt(1+lambda_c^2);
    Delta_c = sigma_c * delta_c;
    Gamma_c = sigma_c^2 * (1-delta_c^2);
    
    zetaO=struct();
    zz=struct();
    zetaO.alpha=alpha; zetaO.beta=beta;
    zetaO.C=struct('mu',u_c,'Gamma',Gamma_c,'Delta',Delta_c);
    zetaO.I1=struct('mu',u_i1, 'Gamma', Gamma_i1, 'Delta', Delta_i1);
    zetaO.I2=struct('mu',u_i2, 'Gamma', Gamma_i2, 'Delta', Delta_i2);
    zetaO.I3=struct('mu',u_i3, 'Gamma', Gamma_i3, 'Delta', Delta_i3);
    zetaN = zetaO;

%     [Vc, Wc] = trunc_norm_moments(delta_c / sigma_c * (S-u_c), sqrt(1-delta_c^2));
%     [Vi, Wi] = trunc_norm_moments(delta_i / sigma_i * (S-u_i), sqrt(1-delta_i^2));
% 
%     Vc1 = Vc(1,:);
%     Vc2 = Vc(2,:);
%     Wc1 = Wc(1,:);
%     Wc2 = Wc(2,:);
%     Vi1 = Vi(1,:);
%     Vi2 = Vi(2,:);
%     Wi1 = Wi(1,:);
%     Wi2 = Wi(2,:);

    [Vc1, Wc1] = trunc_norm_moments(delta_c / sigma_c * (s1-u_c), sqrt(1-delta_c^2));
    [Vc2, Wc2] = trunc_norm_moments(delta_c / sigma_c * (s2-u_c), sqrt(1-delta_c^2));
    [Vi1, Wi1] = trunc_norm_moments(delta_i1 / sigma_i1 * (s1-u_i1), sqrt(1-delta_i1^2));
    [Vi21, Wi21] = trunc_norm_moments(delta_i1 / sigma_i1 * (s2-u_i1), sqrt(1-delta_i1^2));
    [Vi22, Wi22] = trunc_norm_moments(delta_i2 / sigma_i2 * (s2-u_i2), sqrt(1-delta_i2^2));
    [Vi32, Wi32] = trunc_norm_moments(delta_i2 / sigma_i2 * (s3-u_i2), sqrt(1-delta_i2^2));
    [Vi33, Wi33] = trunc_norm_moments(delta_i3 / sigma_i3 * (s3-u_i3), sqrt(1-delta_i3^2));
   
%     S1 = S(1,:);
%     S2 = S(2,:);

    pc1 = skew_norm_pdf(s1, u_c, sigma_c, lambda_c);
    pi11 = skew_norm_pdf(s1, u_i1, sigma_i1, lambda_i1);
    pc2 = skew_norm_pdf(s2, u_c, sigma_c, lambda_c);
    pi21 = skew_norm_pdf(s2, u_i1, sigma_i1, lambda_i1);
    pi22 = skew_norm_pdf(s2, u_i2, sigma_i2, lambda_i2);
    pi32 = skew_norm_pdf(s3, u_i2, sigma_i2, lambda_i2);
    pi33 = skew_norm_pdf(s3, u_i3, sigma_i3, lambda_i3);
    
    p_total = alpha * pc1 .* pi21 .* pi32 + (1-alpha)*beta * pi11 .* pc2 .* pi32 + (1-alpha)*(1-beta) * pi11 .* pi22 .* pi33;
    
    R1 = alpha * pc1 .* pi21 .* pi32 ./ p_total;
    R2 = (alpha*beta * pc1 .* pi21 .* pi32 + (1-alpha)*beta * pi11 .* pc2 .* pi32) ./ p_total;
    
    sum_R1 = sum(R1);
    sum_R2 = sum(R2);
    
    
    alpha_new = sum_R1 / M;
    lla = func_ll3_3(s1, s2, s3, alpha_new, beta, u_c, sigma_c, lambda_c, u_i1, sigma_i1, lambda_i1, u_i2, sigma_i2, lambda_i2, u_i3, sigma_i3, lambda_i3);
    if (lla-prev_ll<0)
        disp('hi')
    end
    zetaN.alpha = alpha_new;
    beta_new = sum_R2 / M;
    llb = func_ll3_3(s1, s2, s3, alpha_new, beta_new, u_c, sigma_c, lambda_c, u_i1, sigma_i1, lambda_i1, u_i2, sigma_i2, lambda_i2, u_i3, sigma_i3, lambda_i3);
    if (llb-lla<0)
        disp('hi')
    end
    zetaN.beta = beta_new;
    u_c_new = (sum( R1 .* (s1 - Vc1 * Delta_c) ) + sum( ((1 - R1).*R2) .* (s2 - Vc2 * Delta_c) )) / (sum_R1 + sum((1 - R1).*R2));
    ll1 = func_ll3_3(s1, s2, s3, alpha_new, beta_new, u_c_new, sigma_c, lambda_c, u_i1, sigma_i1, lambda_i1, u_i2, sigma_i2, lambda_i2, u_i3, sigma_i3, lambda_i3);
    if (ll1-llb<0)
        disp('hi')
    end
    zetaN.C.mu = u_c_new;
    u_i1_new = (sum( R1 .* (s2 - Vi1 * Delta_i1) ) + sum( (1 - R1) .* (s1 - Vi1 * Delta_i1) )) / M;
    ll2 = func_ll3_3(s1, s2, s3, alpha_new, beta_new, u_c_new, sigma_c, lambda_c, u_i1_new, sigma_i1, lambda_i1, u_i2, sigma_i2, lambda_i2, u_i3, sigma_i3, lambda_i3);
    if (ll2-ll1<0)
        disp('hi')
    end
    zetaN.I1.mu = u_i1_new;
    u_i2_new = (sum( (R1 + (1 - R1).*R2) .* (s3 - Vi32 * Delta_i2) ) + sum( ((1 - R1).*(1 - R2)) .* (s2 - Vi22 * Delta_i2) )) / M;
     ll3 = func_ll3_3(s1, s2, s3, alpha_new, beta_new, u_c_new, sigma_c, lambda_c, u_i1_new, sigma_i1, lambda_i1, u_i2_new, sigma_i2, lambda_i2, u_i3, sigma_i3, lambda_i3);
    if (ll3-ll2<0)
        disp('hi')
    end
     zetaN.I2.mu = u_i2_new;
    u_i3_new = (sum( ((1 - R1).*(1 - R2)) .* (s3 - Vi33 * Delta_i3) )) / sum((1 - R1).*(1 - R2));
     ll4 = func_ll3_3(s1, s2, s3, alpha_new, beta_new, u_c_new, sigma_c, lambda_c, u_i1_new, sigma_i1, lambda_i1, u_i2_new, sigma_i2, lambda_i2, u_i3_new, sigma_i3, lambda_i3);
    if (ll4-ll3<0)
        disp('hi')
    end
     zetaN.I3.mu = u_i3_new;
%     ll = func_ll2(S, alpha, u_c_new, Delta_c, Gamma_c, u_i_new, Delta_i, Gamma_i)
%     pc1 = skew_norm_pdf(s1, u_c_new, sigma_c, lambda_c);
%     pi11 = skew_norm_pdf(s1, u_i1_new, sigma_i1, lambda_i1);
%     pc2 = skew_norm_pdf(s2, u_c_new, sigma_c, lambda_c);
%     pi21 = skew_norm_pdf(s2, u_i1_new, sigma_i1, lambda_i1);
%     pi22 = skew_norm_pdf(s2, u_i2_new, sigma_i2, lambda_i2);
%     pi32 = skew_norm_pdf(s3, u_i2_new, sigma_i2, lambda_i2);
%     pi33 = skew_norm_pdf(s3, u_i3_new, sigma_i3, lambda_i3);
%     
%     p_total = alpha_new * pc1 .* pi21 .* pi32 + (1-alpha_new)*beta_new * pi11 .* pc2 .* pi32 + (1-alpha_new)*(1-beta_new) * pi11 .* pi22 .* pi33;
%     
%     R1 = alpha_new * pc1 .* pi21 .* pi32 ./ p_total;
%     R2 = (alpha_new*beta_new * pc1 .* pi21 .* pi32 + (1-alpha_new)*beta_new * pi11 .* pc2 .* pi32) ./ p_total;
%     
%     sum_R1 = sum(R1);
%     sum_R2 = sum(R2);
% [Vc1, Wc1] = trunc_norm_moments(delta_c / sigma_c * (s1-u_c_new), sqrt(1-delta_c^2));
%     [Vc2, Wc2] = trunc_norm_moments(delta_c / sigma_c * (s2-u_c_new), sqrt(1-delta_c^2));
%     [Vi1, Wi1] = trunc_norm_moments(delta_i1 / sigma_i1 * (s1-u_i1_new), sqrt(1-delta_i1^2));
%     [Vi21, Wi21] = trunc_norm_moments(delta_i1 / sigma_i1 * (s2-u_i1_new), sqrt(1-delta_i1^2));
%     [Vi22, Wi22] = trunc_norm_moments(delta_i2 / sigma_i2 * (s2-u_i2_new), sqrt(1-delta_i2^2));
%     [Vi32, Wi32] = trunc_norm_moments(delta_i2 / sigma_i2 * (s3-u_i2_new), sqrt(1-delta_i2^2));
%     [Vi33, Wi33] = trunc_norm_moments(delta_i3 / sigma_i3 * (s3-u_i3_new), sqrt(1-delta_i3^2));
   
    
    Delta_c_new = (sum( R1 .* Vc1 .* (s1 - u_c_new) ) + sum( (1 - R1).*R2 .* Vc2 .* (s2 - u_c_new) )) / (sum(R1 .* Wc1) + sum((1 - R1).*R2 .* Wc2));
     zetaN.C.Delta = Delta_c_new;
    lambda_c_new = sign(Delta_c_new) * sqrt(Delta_c_new^2 / Gamma_c);
    %delta_c_new = sign(Delta_c_new) * sqrt(Delta_c_new^2 / sigma_c^2);
    %lambda_c_new = sign(Delta_c_new)*sqrt(delta_c_new^2/(1-delta_c_new^2));
    sigma_c1=sqrt(Gamma_c + Delta_c_new^2);
     ll5 = func_ll3_3(s1, s2, s3, alpha_new, beta_new, u_c_new, sigma_c1, lambda_c_new, u_i1_new, sigma_i1, lambda_i1, u_i2_new, sigma_i2, lambda_i2, u_i3_new, sigma_i3, lambda_i3);
    if (ll5-ll4<0)
        disp('hi')
    end
    Delta_i1_new = (sum( R1 .* Vi21 .* (s2 - u_i1_new) ) + sum( (1 - R1) .* Vi1 .* (s1 - u_i1_new) )) / (sum(R1 .* Wi21) + sum((1 - R1) .* Wi1));
    zetaN.I1.Delta = Delta_i1_new;
    lambda_i1_new = sign(Delta_i1_new) * sqrt(Delta_i1_new^2 / Gamma_i1);
    sigma_i1_1=sqrt(Gamma_i1 + Delta_i1_new^2);
     ll6 = func_ll3_3(s1, s2, s3, alpha_new, beta_new, u_c_new, sigma_c1, lambda_c_new, u_i1_new, sigma_i1_1, lambda_i1_new, u_i2_new, sigma_i2, lambda_i2, u_i3_new, sigma_i3, lambda_i3);
    if (ll6-ll5<0)
        disp('hi')
    end
    Delta_i2_new = (sum( (R1 + (1 - R1).*R2) .* Vi32 .* (s3 - u_i2_new) ) + sum( (1 - R1).*(1 - R2) .* Vi22 .* (s2 - u_i2_new) )) / (sum((R1 + (1 - R1).*R2) .* Wi32) + sum((1 - R1).*(1 - R2) .* Wi22));
    zetaN.I2.Delta = Delta_i2_new;
     lambda_i2_new = sign(Delta_i2_new) * sqrt(Delta_i2_new^2 / Gamma_i2);
     sigma_i2_1=sqrt(Gamma_i2 + Delta_i2_new^2);
     ll7 = func_ll3_3(s1, s2, s3, alpha_new, beta_new, u_c_new, sigma_c, lambda_c_new, u_i1_new, sigma_i1_1, lambda_i1_new, u_i2_new, sigma_i2_1, lambda_i2_new, u_i3_new, sigma_i3, lambda_i3);
    if (ll7-ll6<0)
        disp('hi')
    end
    Delta_i3_new = (sum( ((1 - R1).*(1 - R2)) .* Vi33 .* (s3 - u_i3_new) )) / (sum( ((1 - R1).*(1 - R2)) .* Wi33));
    zetaN.I3.Delta = Delta_i3_new;
     lambda_i3_new = sign(Delta_i3_new) * sqrt(Delta_i3_new^2 / Gamma_i3);
     sigma_i3_1=sqrt(Gamma_i3 + Delta_i3_new^2);
     ll8 = func_ll3_3(s1, s2, s3, alpha_new, beta_new, u_c_new, sigma_c1, lambda_c_new, u_i1_new, sigma_i1_1, lambda_i1_new, u_i2_new, sigma_i2_1, lambda_i2_new, u_i3_new, sigma_i3_1, lambda_i3_new);
    if (ll8-ll7<0)
        disp('hi')
    end
%     ll = func_ll2(S, alpha, u_c_new, Delta_c_new, Gamma_c, u_i_new, Delta_i, Gamma_i)
    
%     ll = func_ll2(S, alpha, u_c_new, Delta_c_new, Gamma_c, u_i_new, Delta_i_new, Gamma_i)

    Gamma_c_new = ( sum( R1 .* ((s1 - u_c_new).^2 - 2 * Vc1 .* (s1 - u_c_new) * Delta_c_new + Wc1 * Delta_c_new^2) ) ...
        + sum( (1 - R1).*R2 .* ((s2 - u_c_new).^2 - 2 * Vc2 .* (s2 - u_c_new) * Delta_c_new + Wc2 * Delta_c_new^2) ) ) ...
        / (sum_R1 + sum((1 - R1).*R2));
    zetaN.C.Gamma = Gamma_c_new;
     sigma_c_new = sqrt(Gamma_c_new + Delta_c_new^2);
     lambda_c_new = sign(Delta_c_new) * sqrt(Delta_c_new^2 / Gamma_c_new);
      ll9 = func_ll3_3(s1, s2, s3, alpha_new, beta_new, u_c_new, sigma_c_new, lambda_c_new, u_i1_new, sigma_i1, lambda_i1_new, u_i2_new, sigma_i2, lambda_i2_new, u_i3_new, sigma_i3, lambda_i3_new);
    if (ll9-ll8<0)
        disp('hi')
    end
    Gamma_i1_new = ( sum( R1 .* ((s2 - u_i1_new).^2 - 2 * Vi21 .* (s2 - u_i1_new) * Delta_i1_new + Wi21 * Delta_i1_new^2) ) ...
        + sum( (1 - R1) .* ((s1 - u_i1_new).^2 - 2 * Vi1 .* (s1 - u_i1_new) * Delta_i1_new + Wi1 * Delta_i1_new^2) ) ) ...
        / M;
    zetaN.I1.Gamma = Gamma_i1_new;
    sigma_i1_new = sqrt(Gamma_i1_new + Delta_i1_new^2);
    lambda_i1_new = sign(Delta_i1_new) * sqrt(Delta_i1_new^2 / Gamma_i1_new);
      ll10 = func_ll3_3(s1, s2, s3, alpha_new, beta_new, u_c_new, sigma_c_new, lambda_c_new, u_i1_new, sigma_i1_new, lambda_i1_new, u_i2_new, sigma_i2, lambda_i2_new, u_i3_new, sigma_i3, lambda_i3_new);
    if (ll10-ll9<0)
        disp('hi')
    end
    Gamma_i2_new = ( sum( (R1 + (1 - R1).*R2) .* ((s3 - u_i2_new).^2 - 2 * Vi32 .* (s3 - u_i2_new) * Delta_i2_new + Wi32 * Delta_i2_new^2) )...
        + sum( (1 - R1).*(1 - R2) .* ((s2 - u_i2_new).^2 - 2 * Vi22 .* (s2 - u_i2_new) * Delta_i2_new + Wi22 * Delta_i2_new^2) ) ) ...
        / M;
    zetaN.I2.Gamma = Gamma_i2_new;
     sigma_i2_new = sqrt(Gamma_i2_new + Delta_i2_new^2);
     lambda_i2_new = sign(Delta_i2_new) * sqrt(Delta_i2_new^2 / Gamma_i2_new);
      ll11 = func_ll3_3(s1, s2, s3, alpha_new, beta_new, u_c_new, sigma_c_new, lambda_c_new, u_i1_new, sigma_i1_new, lambda_i1_new, u_i2_new, sigma_i2_new, lambda_i2_new, u_i3_new, sigma_i3, lambda_i3_new);
    if (ll11-ll10<0)
        disp('hi')
    end
    Gamma_i3_new = ( sum( (1 - R1).*(1 - R2) .* ((s3 - u_i3_new).^2 - 2 * Vi33 .* (s3 - u_i3_new) * Delta_i3_new + Wi33 * Delta_i3_new^2) ) ) ...
        / sum((1 - R1).*(1 - R2));
    zetaN.I3.Gamma = Gamma_i3_new;
     sigma_i3_new = sqrt(Gamma_i3_new + Delta_i3_new^2);
     lambda_i3_new = sign(Delta_i3_new) * sqrt(Delta_i3_new^2 / Gamma_i3_new);
      ll12 = func_ll3_3(s1, s2, s3, alpha_new, beta_new, u_c_new, sigma_c_new, lambda_c_new, u_i1_new, sigma_i1_new, lambda_i1_new, u_i2_new, sigma_i2_new, lambda_i2_new, u_i3_new, sigma_i3_new, lambda_i3_new);
    if (ll12-ll11<0)
        disp('hi')
    end
%     ll = func_ll2(S, alpha, u_c_new, Delta_c_new, Gamma_c_new, u_i_new, Delta_i_new, Gamma_i)

%     ll = func_ll2(S, alpha, u_c_new, Delta_c_new, Gamma_c_new, u_i_new, Delta_i_new, Gamma_i_new)

    alpha = alpha_new;
    beta = beta_new;
    u_c = u_c_new;
    u_i1 = u_i1_new;
    u_i2 = u_i2_new;
    u_i3 = u_i3_new;
    Delta_c = Delta_c_new;
    Delta_i1 = Delta_i1_new;
    Delta_i2 = Delta_i2_new;
    Delta_i3 = Delta_i3_new;
    Gamma_c = Gamma_c_new;
    Gamma_i1 = Gamma_i1_new;
    Gamma_i2 = Gamma_i2_new;
    Gamma_i3 = Gamma_i3_new;
    lambda_c = sign(Delta_c) * sqrt(Delta_c^2 / Gamma_c);
    lambda_i1 = sign(Delta_i1) * sqrt(Delta_i1^2 / Gamma_i1);
    lambda_i2 = sign(Delta_i2) * sqrt(Delta_i2^2 / Gamma_i2);
    lambda_i3 = sign(Delta_i3) * sqrt(Delta_i3^2 / Gamma_i3);
    sigma_c = sqrt(Gamma_c + Delta_c^2);
    sigma_i1 = sqrt(Gamma_i1 + Delta_i1^2);
    sigma_i2 = sqrt(Gamma_i2 + Delta_i2^2);
    sigma_i3 = sqrt(Gamma_i3 + Delta_i3^2);
    
%     break
%     disp(ll - prev_ll);
%     disp(alpha);
%     sigma_c
%     sigma_i
%     u_c
%     u_i
%     lambda_c
%     lambda_i
end

end
