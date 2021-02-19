function u = searchUB(xx,x)
        l=1;
        u=length(xx);
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
 
  