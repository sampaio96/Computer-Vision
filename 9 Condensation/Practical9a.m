function r = Practical9a
% Overall, you should learn how to use the Condensation algorithm in this
% practical. 
%
% Part a: Factored Sampling
% This part is mostly illustrative, so you just have a few TODO's.
% Condensation is actually Factored Sampling, but applied iteratively to
% a sequence of observations, and incorporating a motion model. So here,
% time will essentially stand still: given a SINGLE set of "observations",
% can you estimate the posterior probabilities? 
% Observations will be simulated: The real 2D distribution we would like 
% to estimate is the red channel of an abstract image (provided). But
% pretend you can't look at the whole image, and can only take measurements
% here and there. 
% Below, you have the code for factored sampling, but note the comment
% "Loop from here". Looping will only be needed in Part b because we will have a
% changing state, so each loop will advance from t to t+1. Here, you can
% abuse factored sampling a little, and loop "in place". If time were
% advancing and you had a real motion model, that would be Condensation.
% Observe: when you do this factored REsampling, more of the particles 
% should be landing near the peaks in the distribution.


img = double(imread( 'abstract.png' ));
MeasurementsComprehensive = img(:,:,1);
imagesc(MeasurementsComprehensive/255)
colormap(gray);
[imgHeight, imgWidth, colors] = size( img );

numParticles = 50;
weight_of_samples = ones(numParticles,1);

% TO DO: normalize the weights (may be trivial this time)
weight_of_samples = weight_of_samples./sum(weight_of_samples);

% Initialize which samples from "last time" we want to propagate: all of
% them!:
samples_to_propagate = [1:numParticles]';



numDims_w = 2;
% Here we randomly initialize some particles throughout the space of w.
% The positions of such particles are quite close to the known initial
% position:
particles_old = rand( numParticles, numDims_w);
particles_old(:,1) = particles_old(:,1) * imgHeight;
particles_old(:,2) = particles_old(:,2) * imgWidth;

hImg = figure;
hSamps = figure;


iTime = 0;
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% After going through the following code once, you can
%       "Loop from here", to the bottom, to see what factored REsampling
%       would look like.
while iTime < 4
iTime = iTime + 1;

% TO DO: compute the cumulative sum of the weights. You could refer to
% the MATLAB built-in function 'cumsum'.
cum_hist_of_weights = cumsum(weight_of_samples);


% ==============================================================
% Resample the old distribution (at time t-1), favoring samples that had
% a higher posterior probability.
% ==============================================================
samples_to_propagate = zeros(numParticles,1);

% Pick random thresholds (numParticles X 1 matrix) in the cumulative
% probability's range [0,1]: 
some_threshes = rand(numParticles,1);


% For each random threshold, find which sample in the ordered set is the
% first one to push the cumulative probability above that threshold. 
% E.g. if the cumulative histogram goes from 0.23 to 0.26 between the 17th
% and 18th samples in the old distribution, and the threshold is 0.234,
% then we'll want to propagate the 18th sample's w (i.e. particle #18).
for (sampNum = 1:numParticles) 
    thresh = some_threshes(sampNum);
    for (index = 1:numParticles)
        if( cum_hist_of_weights(index) > thresh )
            break;
        end;
    end;
    samples_to_propagate(sampNum) = index;
end;
% Note: it's ok if some of the old particles get picked repeatedly, while
% others don't get picked at all.


% =================================================
% Visualize
% =================================================
set(0,'CurrentFigure',hSamps)
set(gcf,'Position',[23 125 640 480]);
set(gcf,'Color',[1 1 1]);
set(gca,'Box','Off');
title('Cumulative histogram of probabilities for sorted list of particles');
plot(zeros(numParticles,1), some_threshes,'b.','MarkerSize',15);
hold on;
plot([1:numParticles], cum_hist_of_weights, 'ro-', 'MarkerSize',3);
which_sample_ids = unique(samples_to_propagate);
how_many_of_each = histc(samples_to_propagate, unique(samples_to_propagate));
for( k=1:size(which_sample_ids,1) )
    plot( which_sample_ids(k), 0, 'bo-', 'MarkerSize', 3 * how_many_of_each(k) )
end;
xlabel(sprintf( 'Indeces of all available samples, with larger blue circles for frequently re-sampled particles\n(Iteration %d)', iTime));
ylabel('Cumulative probability');
hold off
xlim([0 numParticles]);ylim([0 1]);
% =================================================
% =================================================


% Predict where the particles we plucked from the old distribution of 
% state-space will go in the next time-step. This means we have to apply 
% the motion model to each old sample.

particles_new = zeros( size(particles_old) );
for (particleNum = 1:numParticles)
    % TO DO: Incorporate some noise, e.g. Gaussian noise with std 10,
    % into the current location (particles_old), to give a Brownian
    % motion model.
    noise = randn(1,2)*10;
    particles_new(particleNum,:) = particles_old( samples_to_propagate(particleNum),: ) + noise;
end;


set(0,'CurrentFigure',hImg)
imagesc(MeasurementsComprehensive/255)
colormap(gray);
set(gcf,'Position',[23 125 640 480]);
set(gcf,'Color',[1 1 1]);
title(sprintf( 'Particles projected to measurement-space\n(Iteration %d)', iTime));
hold on
plot(particles_new(:,2), particles_new(:,1), 'rx')
hold off


% From here we incorporate the sensor measurement for the new state (time t):
% The new particles, accompanied with predicted locations in world state-space
% for time t, are missing their weights: how well does each particle
% explain the observation x_t?
for (particleNum = 1:numParticles)

    % Convert the particle from state-space w to measurement-space x:
    % Note: It is trivial in this case because both are in 2D space of image
    % coordinates
    
    % Within the loop, we evaluate the likelihood of each particle:
    particle = particles_new(particleNum,:);
    % Check that the predicted location is a place we can really evaluate
    % the likelihood.
    inFrame = particle(1) >= 1.0   &&  particle(1) <= imgHeight && ...
            particle(2) >= 1.0   &&  particle(2) <= imgWidth;
    if( inFrame )
        weight_of_samples(particleNum) = ...
            interp2( MeasurementsComprehensive, ...
            particles_new(particleNum,2), particles_new(particleNum,1) );
    else
        weight_of_samples(particleNum) = 0.0;
    end;

end;

% TO DO: normalize the weights 
weight_of_samples = weight_of_samples./sum(weight_of_samples);


% Now we're done updating the state for time t. 
% For Condensation, just clean up and prepare for the next round of 
% predictions and measurements:
particles_old = particles_new;
clear particles_new;
end
