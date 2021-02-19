function [diff_ss,diff_mu] = calc_diff_kernels(s1, s2, J)
%CALC_KERNELS Summary of this function goes here
%   Detailed explanation goes here
diff = s1-s2; 
diff_ss=nan(length(s2),1);
diff_mu=nan(length(s2),1);
[s2_sorted, s_ind] = sort(s2);
diff_sorted = diff(s_ind);
for i = 1:length(s1)
    ix=find(s2_sorted<s2_sorted(i)+J & s2_sorted>s2_sorted(i)-J) ;
    IX1=[IX1,ix];
    diff_ss(i)= var(diff(ix)); 
    diff_mu(i)= mean(diff(ix)); 
end
        
% for i = 1:length(s1)
%    ix=find(s2<s2(i)+J & s2>s2(i)-J) ;
%    IX1=[IX1,ix];
%    diff_ss(i)= var(diff(ix)); 
%    diff_mu(i)= mean(diff(ix)); 
% end
end

