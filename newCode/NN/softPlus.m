classdef softPlus < nnet.layer.Layer
    % Example custom PReLU layer.

    properties (Learnable)
        % Layer learnable parameters
            
        % Scaling coefficient
    
    end
    
    methods
        function layer = softPlus(numChannels, name) 
            % layer = preluLayer(numChannels, name) creates a PReLU layer
            % with numChannels channels and specifies the layer name.

            % Set layer name.
            layer.Name = name;

            % Set layer description.
            layer.Description = "SoftPlus with " + numChannels + " channels";
        
            % Initialize scaling coefficient.
           
        end
        
        function Z = predict(layer, X)
            % Z = predict(layer, X) forwards the input data X through the
            % layer and outputs the result Z.
            if length(X)>100
                disp('hi');
            end
           % Z = 0.1+log(1+exp(X));
            Z = max(0.1,X);
        end
    end
end