function p = cdfHN(thetaD2,thetaI,ub)
ss=10000;
z = thetaI.mu + thetaI.Delta*abs(randn(ss,1)) + sqrt(thetaI.Gamma)*randn(ss,1);
sigma=sigmaVec(thetaD2,z);
s2 = z - sigma.*abs(randn(ss,1));
p=empiricalCdf(s2,ub);
% xx=sort(xx);
% p=nan(length(ub),1);
% for i=1:length(ub)
%     ubix = searchUB(xx,ub(i));
%     p(i)=ubix/n;
% end
% 
%     function u = searchUB(xx,x)
%         l=1;
%         u=n;
%         while true
%             ix=floor((l+u)/2);
%             if ix==l||ix==u
%                 break;
%             end
%             if xx(ix) < x
%                 l=ix;
%             else
%                 u=ix;
%             end
%         end
%      end
end

