function smp = TNSample(mu,sigSq)
low= -mu./sqrt(sigSq);

n=10000;
xx = randn(n,1);
xx=sort(xx);
smp=nan(length(mu),1);
for i=1:length(mu)
    ubix = searchUB(xx,low(i));
    %lbix = searchLB(xx,low(i));
    ixx=[];
    if xx(ubix)>= low(i)
        ixx = [ixx,ubix:n];
    end
    if ~isempty(ixx)
        ind = ixx(randi(length(ixx)));
        smp(i)=xx(ind);
    else
        
        smp(i)=low(i);
    end
end

smp = mu +  sqrt(sigSq).* smp;

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
 
    function l = searchLB(xx,x)
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

 