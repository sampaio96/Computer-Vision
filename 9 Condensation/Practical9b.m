function r=Practical9b
% Part b: Condensation
%
% Now, we will track a given shape (template) as it moves in a sequence of
% frames. Our shape/appearance model is trivial: just a template, given.
% We will only explore simple motion models. Such models may
% contain a variable that switches between cars and pedestrians, or traffic
% in different directions, or accelerating etc.
% 
% Note: depending on the motion model, the state w will have different
% numbers of dimensions (not just 2, though that's the default). 
%
% Most of this algorithm should be copied from Part a. An important
% difference is that now we'll actually compute the likelihood using the
% (provided) MeasurePatchSimilarityHere() function.
%
%
% Complete remaining TO DO parts you didn't answer in Part a, run the code 
% over the sequence, and observe the results. 
%
% Then TO DO:
% - Try varying the number of particles: 2000, 500, 100,...
% - Change the state to have 2 more degrees of freedom: velX and velY.
%   (This will require you to change how state-predictions are made, and 
%   how they are converted to measurement-space)
% - Visualize the top-scoring particles (more than just 1!)


% Load template and starting position ('pos'), which come from frame 871.
load Template

% Load images from the rest of the sequence
% Warning: this is VERY memory-wasteful! Would normally load-as-we-go.
Imgs = cell(22,1);
iFrame = 0;
for (frameNum = 872:894)
    sImName = sprintf( 'HillsRdSkipFrames_%.7d.png', frameNum );
    iFrame = iFrame + 1;
    Imgs{iFrame} = imread( sImName );
end;
imgWidth = size(Imgs{1}, 2);
imgHeight = size(Imgs{1}, 1);


numParticles = 300;
weight_of_samples = ones(numParticles,1);

% TO DO: normalize the weights (may be trivial this time)
weight_of_samples = weight_of_samples./sum(weight_of_samples);

% Initialize which samples from "last time" we want to propagate: all of
% them!:
samples_to_propagate = [1:numParticles]';


% ============================
% NOT A TO DO: You don't need to change the code below, but eventually you may
% want to vary the number of Dims.
numDims_w = 4;
% Here we randomly initialize some particles throughout the space of w:
particles_old = rand( numParticles, numDims_w);
particles_old(:,1) = particles_old(:,1) * imgHeight;
particles_old(:,2) = particles_old(:,2) * imgWidth;
% ============================


hImg = figure;
hSamps = figure;

% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% Was our "Loop from here" point in Part a.

for( iTime = 1:22 )

    % TO DO: compute the cumulative sume of the weights. You could refer to
    % the MATLAB built-in function 'cumsum'.
    cum_hist_of_weights = cumsum(weight_of_samples);


    % ==============================================================
    % Resample the old distribution at time t-1, and select samples, favoring
    % those that had a higher posterior probability.
    % ==============================================================
    samples_to_propagate = zeros(numParticles,1);
    
    % Pick random thresholds (numParticles X 1 matrix) in the cumulative
    % probability's range [0,1]:
    some_threshes = rand(numParticles,1);


    % For each random threshold, find which sample in the ordered set is
    % the first one to push the cumulative probability above that
    % threshold, e.g. if the cumulative histogram goes from 0.23 to 0.26
    % between the 17th and 18th samples in the old distribution, and the
    % threshold is 0.234, then we'll want to propagate the 18th sample's w
    % (i.e. particle #18).
    
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
    which_sample_ids = which_sample_ids(find(histc(samples_to_propagate, round(unique(samples_to_propagate))*0.1)));
    for( k=1:size(which_sample_ids,1) )
        plot( which_sample_ids(k), 0, 'bo-', 'MarkerSize', 3 * how_many_of_each(k) )
    end;
    xlabel(sprintf( 'Indeces of all available samples, with larger blue circles for frequently re-sampled particles\n(Iteration %d)', iTime));
    ylabel('Cumulative probability');
    hold off
    xlim([0 numParticles]);ylim([0 1]);
    drawnow
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
        particles_new(particleNum,3:4) = particles_old( samples_to_propagate(particleNum),3:4 ) + randn(1,2);
        particles_new(particleNum,1:2) = particles_old( samples_to_propagate(particleNum),1:2 ) + particles_old( samples_to_propagate(particleNum),3:4 ) + noise;
        particles_new(particleNum,1:2) = round(  particles_new(particleNum,1:2)  ); % Round the particles_new to simplify Likelihood evaluation.
    end;
    % TO DO: Not initially, but change the motion model above to have
    % different degrees of freedom, and optionally completely different
    % motion models.

    
    Im2 = double( Imgs{iTime} );
    
    set(0,'CurrentFigure',hImg)
    imagesc(Im2/255)
    colormap(gray);
    set(gcf,'Position',[23 125 640 480]);
    set(gcf,'Color',[1 1 1]);
    title(sprintf( 'Particles projected to measurement-space\n(Time %d)', iTime));
    hold on
    plot(particles_new(which_sample_ids,2), particles_new(which_sample_ids,1), 'rx')
    hold off
    drawnow

    % Optional code to save out figure:
    %     pngFileName = sprintf( '%s_%.5d.png', 'myOutput', iTime );
    %     saveas(gcf, pngFileName, 'png');

    
    % From here we incorporate the data for the new state (time t):
    % The new particles accompanying predicted locations in state-space
    % for time t, are missing their weights: how well does each particle
    % explain the observations x_t?
    for (particleNum = 1:numParticles)

        % Convert the particle from state-space w to measurement-space x:
        % Note: that step is trivial here since both are in 2D space of image
        % coordinates

        % Within the loop, we evaluate the likelihood of each particle:
        particle = particles_new(particleNum,:);
        % Check that the predicted location is a place we can really evaluate
        % the likelihood.
        s = size(pixelsTemplate);
        inFrame = particle(1) >= 1.0   &&  particle(1)+ s(1) <= imgHeight && ...
                particle(2) >= 1.0   &&  particle(2) + s(2) <= imgWidth;
        if( inFrame )
            minX = particles_new(particleNum,2);
            minY = particles_new(particleNum,1);
    
            % A stand-alone demo (ExhaustiveTemplateSearchDemo.m) is
            % provided to illustrate the functionality of our likelihood
            % function:
            % MeasurePatchSimilarityHere(im,template,x,y)
            % im : image to search
            % template : the given template
            % x, y : the upper left corner in the image to place the
            % template for similarity evaluation
            weight_of_samples(particleNum) = ...
                MeasurePatchSimilarityHere( Im2, pixelsTemplate, minY, minX );
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
    particles_old(:,3:4)
    clear particles_new;

end; % End of for loop over each frame in the sequence
