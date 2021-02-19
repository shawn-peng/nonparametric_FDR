function NN = HNlossNN(x,y,p)

X=[ones(length(p),1),x];


function L2 = lossGradient(y,mu,sigmaSq,p)
sigmaSq=par(:,2);
mu=par(:,1);
L1= -0.5 * log (2*pi* sigmaSq) - 0.5 * (y-mu)./sigmaSq;
L2= p' * L1;
end

function l = loss(y,par,p)
sigmaSq=par(:,2);
mu=par(:,1);
L1= -0.5 * log (2*pi* sigmaSq) - 0.5 * (y-mu)./sigmaSq;
l= p' * L1;
end
end