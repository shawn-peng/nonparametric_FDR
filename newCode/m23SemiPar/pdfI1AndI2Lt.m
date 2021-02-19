function p = pdfI1AndI2Lt(zeta,I1,up)
n=10000;
w=nan(length(I1),n);
x=zeta.I2.random(n);
xx=double(x);
[mu,ss]=predictMuAndSS(zeta.D1,x); %% <-- need to use NN for prediction
%p2=snPdf(x,zeta.I2);

% I1_rep = repmat(I1, 1, n);
% up_rep = repmat(up, 1, n);
% xx_rep = repmat(xx', length(I1), 1);
% mu_rep = repmat(mu', length(I1), 1);
% ss_rep = repmat(ss', length(I1), 1);
% diff = I1_rep - xx_rep;
% flags = (up_rep - xx_rep) >= 0;
% w(flags) = tnPdf(diff(flags), mu_rep(flags), sqrt(ss_rep(flags)));
% w = tnPdf(diff, mu_rep, sqrt(ss_rep));
% w = tnPdf((diff-mu_rep)./sqrt(ss_rep), 0, 1);
% w(~flags) = 0;

%% Time previous implementation
for i=1:n
    diff=I1-xx(i);
    jx=diff>=I1-up;
    w(jx,i)= tnPdf(diff(jx),mu(i),sqrt(ss(i)));
    w(~jx,i)=0;
end

%%
p=mean(w,2);
%mn=min(p(p>0));
%p(p==0)=mn;
end
