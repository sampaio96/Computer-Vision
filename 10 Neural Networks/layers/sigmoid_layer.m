% Part 1: sigmoid_layer
% 
% This script contains the class definition for a sigmoid layer. It 
% contains one function only 'forward'.

classdef sigmoid_layer 
    methods
        function [y, obj] = forward(obj, x)
            % Build the forward propagation step for a sigmoid layer.
            y = 1 ./ (1 + exp(-x));
        end
        function [dLdx, obj] = backward(obj, dLdy)
            % Compute the back-propagated gradients of this layer.
            sigmoid = 1 ./ (1 + exp(-x));
            dydx = sigmoid .* (1-sigmoid);
            dLdx = dLdy.*dydx;
            
            % Store gradients to object
            obj.dLdW = [];
        end
    end
end
