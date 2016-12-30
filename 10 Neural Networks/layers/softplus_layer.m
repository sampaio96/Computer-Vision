% Part 1: softplus_layer
% 
% This script contains the class definition for a softplus layer. It
% contains one function only 'forward'.

classdef softplus_layer 
    methods
        function y = forward(obj, x)
            % Build the forward propagation step for a softmax layer.
            % Prevent numerical overflow by subtracting max element of x
            y = log(1 + exp(x));
        end
    end
end