classdef HNLossLayer < nnet.layer.RegressionLayer
   
    
    properties(Constant)
        % Small constant to prevent division by zero.
        Epsilon = 1e-8;
    end
    
    properties
        % Default weighting coefficients 
      
        
    end

    
    methods
        
        function layer = HNLossLayer(name)
            
            % Set layer name.          
            layer.Name = name;
           
            % Set layer description.
            layer.Description = 'HNLoss';
        end
        
        
        function loss = forwardLoss(layer, Y, T)
            sigmaSq=Y(:,:,2,:);
            mu=Y(:,:,1,:);
            %sigmaSq=1;
%             subplot(2,1,1)
%             histogram(sigmaSq, 'Normalization', 'pdf')
%             subplot(2,1,2)
%             histogram(mu, 'Normalization' , 'pdf')
            y=T(:,:,1,:);
            p=T(:,:,2,:);
            L1= 0.5 * log (2*pi* sigmaSq) + 0.5 * ((y-mu).^2)./sigmaSq;
            loss= mean(p .* L1);
            disp(loss);
        end     
    end
end