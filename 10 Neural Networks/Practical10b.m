function Practical10b
% Machine Vision Neural Network tutorial---Part 2
%
% We are going to write a script to run a 7-layer autoencoder. The model is:
%
% input --> encoder --> latent_code --> decoder --> reconstruction
%
% Generate images that look like numbers by generating
% random latent codes from a 500D standard Gaussian and passing these vectors
% through the decoder.
% Find the subspace of the latent
% code, to smoothly interpolate between numbers in
% latent space.

%% Load data and add files
clear
clear 'functions'
addpath('layers')

% Generate data
mnist_images = load('mnist_test');
mnist_images = mnist_images.X;

% Load params
weights = load('weights');

%% Build MLP
% Construct network as an ordered cell array, where each element is a
% layer
encoder = build_encoder(weights);
decoder = build_decoder(weights);

%% Inference

% Plot some input data points.
n_samples = 225;
figure(1)
plot_tiled_array(mnist_images(1:n_samples,:))
title('MNIST examples')


% Generate 225 random codes as a draw from a 500D standard Gaussian.
% The decoder is able to produce convincing images of
% handwritten digits from random Gaussian draws.

latent_code = randn(n_samples,500);%*10; The digits become noisier and less recognizable.

reconstruction = mlp_forward(decoder, latent_code, false);
figure(2)
plot_tiled_array(reconstruction)
title('Decoded randomly generated latent codes')



% Build a linear approximation to the
% data-manifold in the latent space of the autoencoder.
% Now able to smoothly interpolate between digits,
% effectively enforcing a smooth ordering on the data.

% Create sampling grid in a 2D subspace
lim = 3;
image_dims = ceil(sqrt(n_samples));
range = linspace(-lim,lim,image_dims);
[X,Y] = meshgrid(range, range);
sampling_grid = [reshape(X,[n_samples,1]), reshape(Y,[n_samples,1])];

% Generate a random 500D subspace with 2 degrees of freedom:
% Generate two random vectors.

subspace = randn(2, 500); % Not all generated images are digits. Some images in between
% clusters of digits don't look like digits at all.

latent_code = sampling_grid*subspace;

% Reconstruct images from code
reconstruction = mlp_forward(decoder, latent_code, false);
figure(3)
plot_tiled_array(reconstruction)
title('Random subspace')



% Forward pass mnist_images through the encoder to
% have a collection of latent data points.
% Run PCA on the latent codes and keep the first two principal directions.
% A linear manifold is a better approximation to the true data
% manifold.
latent_code = mlp_forward(encoder, mnist_images, false);

% Compute the covariance matrix of the latent codes
latent_cov = cov(latent_code);


% Do PCA on the sampled code:
% Perform the SVD on the covariance matrix and retain the first two
% rows of V (as a column)
[U L Vt] = svd(latent_cov);
firsttworowsV = Vt(:,1:2); 
principal_subspace = firsttworowsV.';


latent_grid = sampling_grid*principal_subspace;
% Reconstruct images from the codes
reconstruction = mlp_forward(decoder, latent_grid, false);

figure(4)
plot_tiled_array(reconstruction)
title('Fitted subspace')


%% Functions
function plot_tiled_array(images)
% Draw images in tiled-array
n_images = size(images, 1);
length = sqrt(size(images,2));
sqrt_nimages = ceil(sqrt(n_images));
tiled_image = zeros(length*sqrt_nimages,length*sqrt_nimages);
for i=1:n_images
    tile = reshape(images(i,:), [length,length])';
    m = length*floor((i-1)/sqrt_nimages)+1;
    n = int32(mod(i-1,sqrt_nimages)*length)+1;
    tiled_image(m:(m+length-1),n:(n+length-1)) = tile;
end
scrsz = get(groot, 'ScreenSize');
imshow(1-tiled_image)
end

% MLP functions
function network = build_encoder(weights)
% Construct a neural network as an ordered cell array. Each element of the
% cell array is a layer. Just make sure that the dimensionality of each
% layer is consistent with its neighbors.

% Declare each layer
affine_e1 = affine_layer(784, 500);
relu_e1 = relu_layer();   
affine_e2 = affine_layer(500, 500);
relu_e2 = relu_layer();
mu_e = affine_layer(500,500);

% Load pretrained weights
affine_e1.W = weights.encoder_W_e1;
affine_e2.W = weights.encoder_W_e2;
mu_e.W = weights.encoder_W_e_mu;

% Build network as ordered cell array
network = {affine_e1, relu_e1, affine_e2, relu_e2, mu_e};
end

function network = build_decoder(weights)
affine_d1 = affine_layer(500, 500);
relu_d1 = relu_layer();   
affine_d2 = affine_layer(500, 500);
relu_d2 = relu_layer();
affine_d3 = affine_layer(500,784);
sigmoid_d = sigmoid_layer();

affine_d1.W = weights.binary_decoder_W_d1;
affine_d2.W = weights.binary_decoder_W_d2;
affine_d3.W = weights.binary_decoder_W_d_mu;

network = {affine_d1, relu_d1, affine_d2, relu_d2, affine_d3, sigmoid_d};
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
end