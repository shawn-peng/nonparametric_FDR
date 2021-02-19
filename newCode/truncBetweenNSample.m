function smp = truncBetweenNSample(mu,sigSq,upper)
up=(upper-mu)./sqrt(sigSq);
low= -mu./sqrt(sigSq);

n=10000;
xx = randn(n,1);
xx=sort(xx);
smp=nan(length(up),1);
for i=1:length(up)
    ubix = searchUB(xx,up(i));
    lbix = searchLB(xx,low(i));
    ixx=[];
    if xx(ubix)>= up(i)
        ixx = [ixx,ubix:n];
    end
    if xx(lbix)<=low(i)
        ixx = [ixx,1:lbix];
    end
    if ~isempty(ixx)
        ind = ixx(randi(length(ixx)));
        smp(i)=xx(ind);
    else
        uplow=[up(i),low(i)];
        smp(i)=uplow(randperm(2,1));
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

 