% Part 1: softmax_layer
% 
% This script contains the class definition for a softmax classifier. It
% contains one function only 'forward'.

classdef softmax_layer 
    methods
        function y = forward(obj, x)
            % Build the forward propagation step for a softmax layer.
            % Prevent numerical overflow by subtracting max element of x
            a = x - max(x,2);
            y = exp(a) ./ repmat(sum(exp(a),2), [1, size(a,2)]);
        end
    end
end
