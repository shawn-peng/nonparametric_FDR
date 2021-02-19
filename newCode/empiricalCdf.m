function p = empiricalCdf(smp,ub)
n=length(smp);
smp=sort(smp);
p=nan(length(ub),1);
for i=1:length(ub)
    ubix = searchUB(smp,ub(i));
    p(i)=ubix/n;
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

