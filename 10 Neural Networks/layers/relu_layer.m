% Part 1: relu_layer
% 
% This script contains the class definition for a rectified linear unit. It
% contains two functions 'forward' and 'backward'. These compute the
% forward- and back-propagation steps respectively.

classdef relu_layer
    % The properties section lists the variables associated with this layer
    % which are stored whenever the forward or backward methods are called.
    properties
        x       % input
        y       % output
        dLdW    % gradient of loss wrt params
    end
    methods
        function [y, obj] = forward(obj, x)
            % Write the forward propagation step for a ReLU_layer
            y = max(zeros,x);
            
            % Save input/output to object properties
            obj.x = x;
            obj.y = y;
        end
        function [dLdx, obj] = backward(obj, dLdy)
            % Compute the back-propagated gradients of this layer.
            % Note that the softmax contains no parameters, so dLdW 
            % is just an empty array
            
            % Compute the gradients wrt the input
            dydx = gt(obj.x,0);
            dLdx = dLdy.*dydx;
            
            % Store gradients to object
            obj.dLdW = [];
        end
    end
end
