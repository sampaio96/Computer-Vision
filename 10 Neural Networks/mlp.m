% Multilayer Perceptron
% 
classdef mlp
    properties
        net % the network
    end
    methods
        
    function [obj] = run(obj, varargin)
    % The structure of the code below is designed to 
    % highlight the modular, reusable nature of neural 
    % network layers. Each layer is defined in a file 
    % with name "<name>_layer.m", for example, the ReLU 
    % layer is in relu_layer.m. Each layer file has at 
    % least 2 'methods'. There is a 'forward' method, 
    % which implements the  forward-propagation step 
    % needed to compute activations, and a 'backward' 
    % method, which implements the back-propagation step. 

    %% Load data and add files 
    cla
    addpath('layers')

    % Generate data
    [X, t] = mlp.generate_data();
    
    % Find the limits of the data and add boundary for plotting
    min_ = min(X,[],1);
    max_ = max(X,[],1);
    diff = max_ - min_;
    min_ = min_ - diff/3;
    max_ = max_ + diff/3;

    % Create a dense grid of testpoints (test_coords), which we shall use to
    % visualize the forward model p(y|x)
    [I,J] = meshgrid(linspace(min_(1), max_(1)), linspace(min_(2), max_(2)));
    test_coordinates = [I(:),J(:)];

    %% Build MLP
    % Parameters for stochastic gradient descent
    minibatch_size = 10;
    initial_learning_rate = 1e-2;

    % Construct the network as an ordered cell array, where each element is a
    % layer
    obj.net = mlp.build_mlp(size(X,2), 250, size(t,2));

    figure(1);
    %% Training loop
    if nargin>1
        step_number = varargin{1};
    else
        step_number = 10000;
    end
    for i=1:step_number
        % Adaptive learning rate satisfying Robbins-Monro conditions
        learning_rate = initial_learning_rate/sqrt(i);

        % Minibatching
        mb = randi(200,1,minibatch_size);
        xmb = X(mb,:);
        tmb = t(mb,:);

        % Complete the forward propagation step
        [logits, obj.net] = mlp.mlp_forward(obj.net, xmb, true);

        % We construct a separate loss layer on the end of the network, which
        % we interchange for the softmax at test-time. By merging the
        % crossentropy and softmax layers at training time, we avoid numerical
        % precision problems during gradient descent.

        % Complete the forward propagration step
        loss_layer = crossentropy_softmax_layer();
        [loss, loss_layer] = forward(loss_layer, logits, tmb);

        % Implement the backward pass
        [dLdy, ~] = backward(loss_layer, 1);
        obj.net = mlp.mlp_backward(obj.net, dLdy);

        % Implement stochastic gradient descent code
        obj.net = mlp.apply_gradient_descent_step(obj.net, learning_rate);

        % Validation
        % Get training accuracy on this minibatch
        [~,indy] = max(logits,[],2);
        [~,indt] = max(tmb,[],2);
        accuracy = mean(indt==indy);

        % Test
        fprintf('[%04i], Loss: %f, Accuracy: %f, Learning rate: %f\n', i, loss, accuracy, learning_rate);
        % Run test data through mlp
        test_logits = mlp.mlp_forward(obj.net, test_coordinates, false);
        test_output = forward(softmax_layer, test_logits);
        % Plot mesh
        pcolor(I,J,reshape(test_output(:,1),size(I)));
        hold on;
        scatter(X(:,1), X(:,2), 10, [t, ones(200,1)]);
        shading flat;
        hold off;
        pause(0.001)
    end
    end
    end
    
    methods(Static)

    function [X, t] = generate_data()
    mean1 = [0,0];
    std1 = 0.2;
    mean2 = [0,1];
    std2 = 0.2;
    % We shall generate two interlocking arcs
    theta1 = 1.2*pi*(rand(100,1)-0.5);
    X1 = repmat(mean1,[100,1]) + [cos(theta1),sin(theta1)] + std1*randn(100,2);
    theta2 = 1.2*pi*(rand(100,1)+0.5);
    X2 = repmat(mean2,[100,1]) + [cos(theta2),sin(theta2)] + std2*randn(100,2);

    % Input X and target t
    X = [X1; X2];
    t = [ones(100,1), zeros(100,1);
        zeros(100,1), ones(100,1)];
    end


    %% Functions
    function network = build_mlp(n_in, n_hid, n_out)
    % Construct a neural network as an ordered cell array. Each element of the
    % cell array is a layer. Just make sure that the dimensionality of each
    % layer is consistent with its neighbors.

    % Declare each layer
    affine1 = affine_layer(n_in, n_hid);
    relu1 = relu_layer();   % ReLU doesn't alter dimensions -> no dim. args
    affine2 = affine_layer(n_hid, n_out);
    

    % Build network as ordered cell array
    network = {affine1, relu1, affine2};
    end

    function [y, net] = mlp_forward(net, x, train)
    % Forward-propagation has 2 modes:
    %   train=true: store the results of the forward pass in each layer
    %   train=false: do not store the results of the forward pass
    % Each layer takes as input the output y from the layer below.

    y = x;
    if train == true
        for j=1:length(net)
            [y, net{j}] = net{j}.forward(y);
        end
    else
        for j=1:length(net)
            [y, ~] = net{j}.forward(y);
        end
    end
    end

    function net = mlp_backward(net, dLdy)
    % Back-propagation: Each layer takes as input the back-propagated 
    % errors/gradients/deltas df from the layer above.
    %
    % CAVEAT: Gradient computation relies not only on dLdy but also x and y 
    % from forward-propagation step, these are stored in-place in each layer, 
    % so you can only run this method after calling 
    % mlp_forward(net, x, train=true).

    % Implement the backward pass
    for j=1:length(net)
        [dLdy, net{end-j+1}] = net{end-j+1}.backward(dLdy);
    end
    end

    function net = apply_gradient_descent_step(net, learning_rate)
    % Gradient descent step: Apply simple stochastic gradient descent step
    for j=1:length(net)
        if ~isempty(net{j}.dLdW)
            % TODO 4: update the weights of the multilayer perceptron
            net{j}.W = net{j}.W - learning_rate*net{j}.dLdW;
        end
    end
    end

    end
end