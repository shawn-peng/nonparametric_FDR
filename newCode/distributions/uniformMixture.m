classdef uniformMixture < mixture 
    methods
        function dens = pdf(this,x)
            if ~isa(x,'binnedDouble')
                dens=this.pdf@mixture(x);
            else
                xx=double(x);
                ixBinned=x.IX;
                dens=nan(length(x),1);
                for i=1:length(this.mixProp)
                    ix=ixBinned{i};
                    if(~isempty(ix))
                        pp=pdf(this.comps{i},xx(ix(1)));
                        dens(ix)=this.mixProp(i)*pp;
                    end
                end
            end
        end
        function xBinned = bin(this,x)
            x=x(:);
            IX=cell(length(this.mixProp),1);
            ix = 1:length(x);
            ix=ix(:);
            for i=1:length(this.mixProp)
                pp=pdf(this.comps{i},x);
                IX{i}=ix(pp>0);
            end
            xBinned=binnedDouble(x,IX);
        end
        
        function mix = fixedCompsFit(this,x,weights)
            if ~isa(x,'binnedDouble')
               mix = this.fixedCompsFit@mixture(x);
            else
                if nargin < 3
                    weights = ones(length(x),1);
                end
                ixBinned=x.IX;
                for i=1:length(this.mixProp)
                    ix=ixBinned{i};
                    if any(ix>length(weights))
                        disp('something is wrong')
                    end
                    this.mixProp(i)=sum(weights(ix));
                end
                this.mixProp=this.mixProp/sum(this.mixProp);
                mix=this;
            end
        end
        
        function c = cdf(this,x)
            ixBinned=x.IX;
            xx=double(x);
            c= zeros(length(x),1);
            ix = 1:length(x);
           
            for i = 1: length(this.mixProp)
                jx=ixBinned{i};
                ix= setdiff(ix,jx);
                c(ix)= c(ix) + this.mixProp(i);
                c(jx) = c(jx) + this.mixProp(i)*cdf(this.comps{i},xx(jx));
            end
        end
        
        function xBinned = random(this,n)
            ncomps=length(this.mixProp);
            x=nan(n,1);
            compSS=floor(n*this.mixProp);
            compSS(ncomps)=compSS(ncomps)+ n -sum(compSS);
            ind=0;
            IX=cell(ncomps,1);
            for i = 1:ncomps
                ix=(ind+1):(ind+compSS(i));
                IX{i}=ix(:);
                x(ix)=random(this.comps{i},compSS(i),1);
                ind=ind+compSS(i);
            end
            xBinned=binnedDouble(x,IX);
        end
        
        function xBinned = ltRandom(this,ub)
            mp=this.mixProp;
            ubb=double(ub);
            mp=mp(:);
            ncomps=length(mp);
            x=nan(length(ub),1);
            IX=cell(ncomps,1);
            for i = 1:ncomps
                low=this.comps{i}.Lower;
                ix=ub.IX{i};
                p = mp(i)*cdf(this.comps{i},ubb(ix));
                for k=1:length(ix)
                    j=ix(k);
                    wt=[mp(1:(i-1));p(k)];
                    if (sum(wt)>0)
                        ii=randsample(i,1,true,wt);
                        if (ii == i)
                            x(j)=rand(1)*(ubb(j)-low)+ low;
                        else
                            x(j)=this.comps{ii}.random(1);
                        end
                    else
                        ii=i;
                        x(j)=ubb(j);
                    end
                        IX{ii}=[IX{ii};j];
                end
            end
            xBinned=binnedDouble(x,IX);
        end
  
   end
end

