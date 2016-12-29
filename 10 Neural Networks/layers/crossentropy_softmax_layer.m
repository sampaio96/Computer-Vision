% Part 1: crossentropy_softmax_layer
% 
% This script contains the class definition for the crossentropy and
% softmax layer combined. This is a common techinque, used to avoid 
% numerical overflow in backpropagation. It contains two functions 
% 'forward' and 'backward'. These compute the forward- and back-propagation 
% steps respectively. You will need to fill out the sections marked 'TODO'.

classdef crossentropy_softmax_layer
    % The properties section lists the variables associated with this layer
    % which are stored whenever the forward or backward methods are called.
    properties
        softmax_output  % Softmax output of network
        target          % Targets
        x               % Input
        y               % Output loss
        dLdW            % Gradient of loss wrt params
    end
    methods
        function [y, obj] = forward(obj, x, target)
            % Compute the forward-propagated activations of this layer. You
            % can do this by computing the softmax followed by the
            % crossentropy.
            
            % subtract max(x) for stability (note there is a MATLAB
            % version, which is defined differently!)
            xx = x - max(x,2);
            obj.softmax_output = softmax(obj, xx);
            
            % TODO 3.2: Compute the crossentropy
            y = -sum(sum(target .* log(obj.softmax_output) + (1-target) .* log(1 - obj.softmax_output))) / size(target(:),1);
            
            % Store the input, output and target to object
            obj.x = x;
            obj.y = y;
            obj.target = target;
        end
        function [dLdx, obj] = backward(obj, dLdy)
            % Compute the back-propagated gradients of this layer.
            % Note that the softmax contains no parameters, so dLdW is
            % empty. Also note that the input gradient dLdy is a scalar.
            
            % TODO 3.3: Compute the gradients wrt the input. (Cheatsheet)
            dLdx = obj.softmax_output-obj.target;
            
            % Store gradients to object
            obj.dLdW = [];
        end
        function y = softmax(obj, xx)
            % TODO 3.1: Compute softmax output
            y = exp(xx) ./ repmat(sum(exp(xx),2), [1, size(xx,2)]);
        end
    end
end
