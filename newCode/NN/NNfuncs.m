function [init, run, predict] = NNfuncs(WIDTH,ACTS,loss)
init=@initialize;
run=@runn;
predict=@predictt;

        function par = initialize()
            for i=1:(length(WIDTH)-1)
                %XXX=randn(WIDTH(i+1),100);
                %[~,C]=kmeans(XXX',WIDTH(i)+1);
                %Ws{i}=sqrt(2/WIDTH(i))*C';
                par.Ws{i}=2*randn(WIDTH(i+1),WIDTH(i)+1);
                par.BNL{i}= randn(WIDTH(i),1);
                par.BNS{i}= randn(WIDTH(i),1);
            end
           % Ws{1}=[1,0];
           %Ws{2}=[1,0];
        end
        
    function Ws= runn(par,x,epoch,eta)
        losses=[];
        %X=[x';ones(1,length(x))];
        Ws=par.Ws;
        BNS=par.BNS;
        BNL=par.BNL;
        for k=1:epoch
            ind=getIX(length(x));
            %Xk=X(:,ind);
            xk=x(ind)';
            %yk=y(ind);
            %pk=p(ind);
            %vk=v(ind);
            
            %IN={Xk};
            F={};
            G={};
            f=xk;
            IN={};
            %BN={};
            n=length(xk);
            for i=1:(length(WIDTH)-1)
                xx=f;
                xx=xx-repmat(mean(xx,2),n);
                S=diag(1./sqrt(var(xx,[],2)));
                xx=S*xx;
                xx = diag(BNS{i})*xx + repmat(BNL{i},n);
                IN{i}=[xx,ones(1,n)];
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
            disp(f);
            losses=[losses,f];
            if k>100
                if losses(k-100)-losses(k)<10^-3
                    eta = max(eta/2,0.001);
                end
            end
            if any(isnan(g),'all')
                disp('nan Produced');
            end
            WW=Ws;
            del={};
            del{i}=g.*G{i};
            WW{i}=outProdSum(del{i},IN{i});
            
            
            for ix=(i:-1:2)
                w=Ws{ix};
                w=w(:,1:(size(Ws{ix},2)-1));
                del{ix-1} = (w'*del{ix}).* G{ix-1};
                WW{ix-1}=outProdSum(del{ix-1},IN{ix-1});
            end
           newWs=Ws;
            for ix=1:(length(WIDTH)-1)
               eta=1;
               while eta~=0
                    newWs{ix}=Ws{ix};
                    newWs{ix}= newWs{ix} - eta*WW{ix};
                    ff=predictt(newWs,x(ind));
                    [ff,~]=loss(ff',ind);
                    if ff > f
                       eta=eta/2;
                    else
                        %disp(ff<f);
                        break; 
                    end
                end
            end
            Ws=newWs;
        end
    end
        
     function yhat = predictt(Ws,x)
        XX=[x';ones(1,length(x))];
        INP{1}=XX;
        for ii=1:(length(WIDTH)-1)
            f=Ws{ii}*INP{ii};
            act=ACTS{ii};
            [f,~]=act(f);
            if ii<(length(WIDTH)-1)
                INP{ii+1}=[f;ones(1,size(f,2))];
            else
                INP{ii+1}=f;
            end
        end
        yhat=f(1,:)';
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
        ix=randi(n,1000,1);
    end
end


