% Machine Vision Neural Network tutorial---Part 1
% 
% You are going to train a small multilayer perceptron
% (MLP) to classify some toy data points in a 2D-space. 
% 
% The multilayer perceptron is implemented in "mlp.m" file..

% Below, the layer functions are tested to ensure that the 
% implementation is correct.

% Make layer files callable.
addpath('layers')
% Fix random seed
rng(1337)

clear;
relu = relu_layer();

% Implement forward pass of ReLU layer
test_x = [-0.1, 0.1; 0.1, -0.1];
[y, relu] = relu.forward(test_x);
test_y = [0.0, 0.1; 0.1, 0.0];
if numel(y)~=numel(test_y) || any(abs(y(:)-test_y(:))>1.0e-6)
    error('TODO 1.1 ReLU forward pass test failed.')
else
    disp('ReLU forward pass test correct.')
end

% Implement backward pass of ReLU layer
test_dldy = [0.1, 0.2; 0.3, 0.4];
[dldx, relu] = relu.backward(test_dldy);
test_dldx = [0.0, 0.2; 0.3, 0.0];
if numel(dldx)~=numel(test_dldx) || any(abs(dldx(:)-test_dldx(:))>1.0e-6)
    error('TODO 1.2 ReLU backward pass test failed.')
else
    disp('ReLU backward pass test correct.')
end

clear;
affine = affine_layer(3, 2);

% Implement forward pass of affine layer
test_x = [-0.1, 0.2, 0.3; 0.1, -0.2, -0.3];
[y, affine] = affine.forward(test_x);
test_y = [-0.15007770, -0.021137288; 0.17007770, 0.041137288];
if numel(y)~=numel(test_y) || any(abs(y(:)-test_y(:))>1.0e-6)
    error('TODO 2.1 Affine forward pass test failed.')
else
    disp('Affine forward pass test correct.')
end

% Implement backward pass of affine layer
test_dldy = [-0.1, -0.2; 0.3, 0.4];
[dldx, affine] = affine.backward(test_dldy);
test_dldx = [0.05000650, 0.08486181, 0.03421172; -0.14523331, -0.21020669, -0.10986740];

if numel(dldx)~=numel(test_dldx) || any(abs(dldx(:)-test_dldx(:))>1.0e-6)
    error('TODO 2.2 Affine backward pass test failed.')
else
    disp('Affine backward pass test correct.')
end

clear;
crossentropy_softmax = crossentropy_softmax_layer();

% Implement "softmax" function of cross-entropy softmax layer 
test_x = [-0.1, 0.3; 0.2, -0.3];
test_target = [0.0, 1.0; 1.0, 0.0];
[y, crossentropy_softmax] = crossentropy_softmax.forward(test_x, test_target);

test_softmax_output = [0.401312339887548, 0.598687660112452; 0.622459331201855, 0.377540668798146];
softmax_output = crossentropy_softmax.softmax_output;

if numel(softmax_output)~=numel(test_softmax_output) || any(abs(softmax_output(:)-test_softmax_output(:))>1.0e-6)
    error('TODO 3.1 Cross-entropy softmax layer''s softmax output test failed.')
else
    disp('Cross-entropy softmax layer''s softmax output test correct.')
end

% Implement crossentropy computation of 
% the forward pass of cross-entropy softmax layer 
test_y = 0.49354611;
if numel(y)~=numel(test_y) || any(abs(y(:)-test_y(:))>1.0e-6)
    error('TODO 3.2 Cross-entropy softmax forward pass test failed.')
else
    disp('Cross-entropy softmax forward pass test correct.')
end

% Implement backward pass of cross entropy softmax layer 
test_dldy = -0.1;
[dldx, crossentropy_softmax] = crossentropy_softmax.backward(test_dldy);
test_dldx = [0.401312339887548, -0.401312339887548; -0.377540668798145, 0.377540668798146];

if numel(dldx)~=numel(test_dldx) || any(abs(dldx(:)-test_dldx(:))>1.0e-6)
    error('TODO 3.3 Cross-entropy softmax backward pass test failed.')
else
    disp('Cross-entropy softmax backward pass test correct.')
end

% Implement "apply_gradient_descent_step" function (backpropagation) 
% of the multilayer perceptron training.
clear;
mlp = mlp();
mlp = mlp.run(5);
test_mlp = load('test_net.mat');

if any(cell2mat(cellfun(@(x,y) isprop(x, 'W') && numel(x.W)~=numel(y.W), test_mlp.net, mlp.net, 'UniformOutput', false))) || ...
        any(cell2mat(cellfun(@(x,y) isprop(x, 'W') && any(x.W(:)-y.W(:)>1.0e-6), mlp.net, test_mlp.net, 'UniformOutput', false)))
    error('TODO 4 MLP apply_gradient_descent_step test failed.')
else
    disp('MLP apply_gradient_descent_step test correct.')
end


clear;
mlp = mlp();
mlp = mlp.run();
