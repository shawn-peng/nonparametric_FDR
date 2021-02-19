function smp = truncSNSample1(theta,lb,ub)
n=10000;
xx = theta.mu + theta.Delta*abs(randn(n,1)) + sqrt(theta.Gamma)*randn(n,1);
xx=sort(xx);
smp=nan(length(lb),1);
for i=1:length(lb)
    lbix = searchLB(xx,lb(i));
    ubix = searchUB(xx,ub(i));
    ixx = randi(ubix-lbix+1);
    ixx = ixx + lbix;
    smp(i)=xx(ixx);
end

    function l = searchLB(xx,x)
        l=1;
        u=n;
        while true
            ix=floor((l+u)/2);
            if ix==l||ix==u
                break;
            end
            if xx(ix) <= x
                l=ix;
            else
                u=ix;
            end
        end
        
    end
    function u = searchUB(xx,x)
        l=1;
        u=n;
        while true
            ix=floor((l+u)/2);
            if ix==l||ix==u
                break;
            end
            if xx(ix) < x
                l=ix;
            else
                u=ix;
            end
        end
     end
end

 