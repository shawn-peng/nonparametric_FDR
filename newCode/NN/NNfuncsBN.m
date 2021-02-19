function [init, run, predict] = NNfuncs(WIDTH,ACTS,loss)
init=@initialize;
run=@runn;
predict=@predictt;
EPS=10^-7;

        function par = initialize()
            for i=1:(length(WIDTH)-1)
                %XXX=randn(WIDTH(i+1),100);
                %[~,C]=kmeans(XXX',WIDTH(i)+1);
                %Ws{i}=sqrt(2/WIDTH(i))*C';
                par.Ws{i}=2*randn(WIDTH(i+1),WIDTH(i)+1);
                par.BNL{i}= randn(WIDTH(i),1);
                par.BNS{i}= abs(randn(WIDTH(i),1));
            end
           % Ws{1}=[1,0];
           %Ws{2}=[1,0];
        end
        
    function [par, badInit]= runn(par,x,epoch,eta,lt)
        losses=[];
        %X=[x';ones(1,length(x))];
        badInit=false;
        for k=1:epoch
            Ws=par.Ws;
            BNS=par.BNS;
            BNL=par.BNL;
            ind=getIX(length(x));
            %Xk=X(:,ind);
            xk=x(ind)';
            n=length(xk);
            %yk=y(ind);
            %pk=p(ind);
            %vk=v(ind);
            
            %IN={Xk};
            F={};
            G={};
            f=xk;
            IN={};
            INN={};
            XX={};
            %BN={};
            n=length(xk);
            for i=1:(length(WIDTH)-1)
                XX{i}=f;
                xx=f;
                M{i}=mean(xx,2);
                par.M=M;
                xx=xx-repmat(M{i},1,n);
                V=var(xx,[],2)+EPS;
                S{i}=1./sqrt(V);
                par.S=S;
                INN{i}=diag(S{i})*xx;
                xx = diag(BNS{i})*INN{i} + repmat(BNL{i},1,n);
                IN{i}=[xx;ones(1,n)];
                f=Ws{i}*IN{i};
                act=ACTS{i};
                [f,g]=act(f);
                F{i}=f;
                G{i}=g;
                if any(isnan(f),'all')
                    disp('nan Produced');
                end
                %     if(i<=2)
                %         subplot(3,2,2*(i-1)+1)
                %         scatter(x,f(1,:))
                %         hold on;
                %         subplot(3,2,2*i)
                %         scatter(x,f(2,:))
                %         hold on;
                %     else
                %         subplot(3,2,5)
                %         scatter(x,f)
                %         hold on;
                %     end
                %     if k==1
                %         subplot(3,2,6)
                %         histogram(f,'Normalization','pdf');
                %     end
                
%                 if i<(length(WIDTH)-1)
%                     IN{i+1}=[f;ones(1,size(f,2))];
%                     %G{i}=[g;zeros(1,size(g,2))];
%                 else
%                     IN{i+1}=f;
%                 end
            end
           
            % subplot(2,2,3)
            % scatter(x,f)
            % hold on;
            [f,g]=loss(f,ind);
            if rem(k,10)==9
                disp(f);
            end
            losses=[losses,f];
%             if k>100
%                 if losses(k-100)-losses(k)<10^-3
%                     eta = max(eta/2,0.001);
%                 end
%             end
            if (k>15 && f>lt) || isnan(f) || isinf(f)
%                 disp('reinitializing NN weights');
                badInit=true;
                break
            end
                
            if any(isnan(g),'all')
                disp('nan Produced');
            end
            WW=Ws;
            del={};
            del{i}=g.*G{i};
            WW{i}=outProdSum(del{i},IN{i});
            
            
            for ix=(i:-1:2)
                xxx=INN{ix};
                w=Ws{ix};
                w=w(:,1:(size(Ws{ix},2)-1));
                W1=w'*del{ix};
                BNLL{ix}=sum(W1,2);
                BNSS{ix}=sum(W1.*xxx,2);
                w=w.*repmat(BNS{ix}',size(w,1),1);
                W2=w'*del{ix};
                tmp = diag(S{ix})*W2;
                Dmu=-(1/n)*sum(tmp,2);
                Dsig= -(1/n)*repmat(sum(tmp.*xxx,2),1,n).*xxx;
                Dx= tmp + Dmu + Dsig;
                del{ix-1} = Dx.* G{ix-1};
                WW{ix-1}=outProdSum(del{ix-1},IN{ix-1});
            end
            xxx=INN{1};
            w=Ws{1};
            w=w(:,1:(size(Ws{1},2)-1));
            W1=w'*del{1};
            BNLL{1}=sum(W1,2);
            BNSS{1}=sum(W1.*xxx,2);
           newWs=Ws;
           newBNS=BNS;
           newBNL=BNL;
           newPar=par;
            for ix=1:(length(WIDTH)-1)
               eta=1;
               while eta~=0
                    newWs{ix}=Ws{ix};
                    newWs{ix}= newWs{ix} - eta*WW{ix};
                    newPar.Ws=newWs;
                    ff=predictt(newPar,x(ind));
                    [ff,~]=loss(ff',ind);
                    if ff > f
                       eta=eta/2;
                    else
                        %disp(ff<f);
                        break; 
                    end
               end
               eta=1;
               while eta~=0
                    newBNL{ix}=BNL{ix};
                    newBNS{ix}=BNS{ix};
                    newBNL{ix}= newBNL{ix} - eta*BNLL{ix};
                    newBNS{ix}= newBNS{ix} - eta*BNSS{ix};
                    newPar.BNS=newBNS;
                    newPar.BNL=newBNL;
                    ff=predictt(newPar,x(ind));
                    [ff,~]=loss(ff',ind);
                    if ff > f
                       eta=eta/2;
                    else
                        %disp(ff<f);
                        break; 
                    end
                end
            end
            par=newPar;
        end
        
    end

    function yhat = predictt(par,x)
        Ws=par.Ws;
        BNL=par.BNL;
        BNS=par.BNS;
        S=par.S;
        M=par.M;
        xk=x';
        %n=length(xk);
        F={};
        G={};
        f=xk;
        IN={};
        INN={};
        XX={};
        %BN={};
        n=length(xk);
        for i=1:(length(WIDTH)-1)
            XX{i}=f;
            xx=f;
            xx=xx-repmat(M{i},1,n);
            INN{i}=diag(S{i})*xx;
           if(size(INN{i},1)~=length(BNS{i}))
               disp('error');
           end
            xx = diag(BNS{i})*INN{i} + repmat(BNL{i},1,n);
            IN{i}=[xx;ones(1,n)];
            f=Ws{i}*IN{i};
            act=ACTS{i};
            [f,g]=act(f);
            F{i}=f;
            G{i}=g;
        end
        yhat=f';
    end

   function ops = outProdSum(v1,v2)
        ops = zeros(size(v1,1),size(v2,1));
        n=size(v1,2);
        for iii=1:n
            ops = ops + v1(:,iii)*v2(:,iii)';
        end
        ops=ops;
   end
    function ix = getIX(n)
        %ix=1:n;
        ix=randi(n,1000,1);
    end
end


