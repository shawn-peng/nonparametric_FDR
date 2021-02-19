function smp = nonTruncNSample(mu,sigSq,x)
low= -mu./sqrt(sigSq);

n=10000;
xx = randn(n,1);
xx=sort(xx);
smp=nan(length(low),1);
for i=1:length(low)
    lbix = searchLB(xx,low(i));
    ixx=[];
    if xx(lbix)<=low(i)
        ixx = [ixx,1:lbix];
    end
   
    if ~isempty(ixx)
        B=binornd(1,length(ixx)/n,1);
        ind=ixx(randi(length(ixx)));
        smp(i)= mu(i) +  sqrt(sigSq(i)).* xx(ind);
        smp(i)=B*smp(i) + (1-B)*x(i);
    else
        smp(i)=x(i);
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

 