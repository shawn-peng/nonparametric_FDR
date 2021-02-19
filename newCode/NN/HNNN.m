function [netMu,netSigmaSq] = HNNN(x,y,p)
% net = network(1,3,[1;1;0],[1; 0; 0],[0 0 0; 1 0 0; 0 1 0],[0 0 1]);
% net.layers{1}.size = 5;
% net.layers{2}.size = 2;
% net.layers{3}.size = 1;
options = trainingOptions('sgdm', ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',0.2, ...
    'LearnRateDropPeriod',5, ...
    'MaxEpochs',20, ...
    'MiniBatchSize',64, ...
    'Plots','training-progress');
netMu.Layers = [ ...
     imageInputLayer([1 1 1])
     fullyConnectedLayer(10)
     %tanhLayer
     reluLayer
     %batchNormalizationLayer
     fullyConnectedLayer(10)
     %tanhLayer
     reluLayer
     %batchNormalizationLayer
     fullyConnectedLayer(10)
     %tanhLayer
     %batchNormalizationLayer
     fullyConnectedLayer(3)
     softPlus(1,'softPlus')
     HNLossLayerMu('HNLoss')
   ];

netSigmaSq.Layers = [ ...
     imageInputLayer([1 1 1])
     fullyConnectedLayer(10)
     %tanhLayer
     reluLayer
      %batchNormalizationLayer
     fullyConnectedLayer(10)
     %tanhLayer
     reluLayer
      %batchNormalizationLayer
     fullyConnectedLayer(10)
     %tanhLayer
      %batchNormalizationLayer
     fullyConnectedLayer(3)
     softPlus(1,'softPlus')
     HNLossLayerSigmaSq('HNLoss')
   ];
XNew = reshape(x', [1,1,size(x,2),size(x,1)]);
muEst=ones(length(x),1);
sigmaSqEst=ones(length(x),1);

for i= 1:2
T=[y,p,sigmaSqEst];
netMu = trainNetwork(XNew, T, netMu.Layers, options);
muEst = predict(netMu,XNew);
muEst=muEst(:,1);
T=[y,p,muEst];
netSigmaSq = trainNetwork(XNew, T, netSigmaSq.Layers, options);
sigmaSqEst = predict(netSigmaSq,XNew);
sigmaSqEst=sigmaSqEst(:,1);
end
end

