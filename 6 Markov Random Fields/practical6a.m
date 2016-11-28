function r= practical6a
%The goal of this practical is to investigate Markov
%random field models.  In practical 6a, we investigate a small 
%5 link model, like the 4-node grid example (12.1.1 in the book).
%
%In 6b, experiment with sampling from probability distributions using 
%MCMC methods (Markov chain Monte Carlo). In practical 6c, we use what
%we have learned to sample from an MRF.


% Note!! Gibbs sampling is illustrated in the book in Figure 10.12, and
% an algorithm outline is provided in the GibbsSampling_Cookbook.pdf, in
% case you want it.


%You should complete the sections of code marked "TO DO".

%close all previous figures
close all;

%make space for output
prW = zeros(1,32);


%run through each of the 32 possible combinations
wIndex = 1;
binaryValues = zeros(32,5);
for(w1 = 0:1)
    for(w2=0:1)
        for(w3=0:1)
            for(w4=0:1)
                for(w5=0:1)
                    %TO DO 
                    %compute unnormalized probability
                    %use the routine getPhi (it's below, also needs completion)
                    %replace this:
                    prW(wIndex) = getPhi(w1,w2)*getPhi(w2,w3)*getPhi(w3,w4)*getPhi(w4,w5);%*getPhi(w5,w1);
                    %store binary values
                    binaryValues(wIndex,:) = [w1 w2 w3 w4 w5];
                    %update index
                    wIndex = wIndex+1;
                end
            end
        end
    end
end

%TO DO
%compute normalization term Z
%replace this
Z = sum(prW);

%normalize probabilities
prW = prW/Z;

%display probabilities
%run through each of the 32 possible combinations
for(cW = 1:32)
    fprintf('Configuration: %d%d%d%d%d, Probability %4.4e\n',binaryValues(cW,1),...
              binaryValues(cW,2),binaryValues(cW,3),binaryValues(cW,4),binaryValues(cW,5),prW(cW));
end;

%TO DO
% Some combinations may appear to have zero probability, when it's really 
% higher probability than zero. Fix the above print statement to make this
% clearer.



% TO DO: fix sampleFromDiscrete() below.
%draw 10 samples from this probability distribution
for (cSample = 1:10)
    w = sampleFromDiscrete(prW);
    fprintf('Sample %d, Configuration: %d%d%d%d%d\n',cSample, binaryValues(w,1), ...
             binaryValues(w,2),binaryValues(w,3),binaryValues(w,4),binaryValues(w,5));    
end;
% Keep going...
%
%
% TO DO: 
% Now change the graphical model and run the above algorithm again: add an
% undirected edge between w5 and w1.

%==========================================================================
%==========================================================================

%return the function phi
function r = getPhi(label1, label2)
%TO DO
%fill in this function with the potentials as given in the smoothing example in the
%notes.
%replace this:
if label1 == label2
    r = 1;
else
    r = 0.1;
end

%==========================================================================
%==========================================================================

%TO DO (see inside function)
%draws a random sample from a discrete probability distribution using a
%rejection sampling method
function r = sampleFromDiscrete(probDist);

nIndex = length(probDist);
while(1)
    %choose random index
    r=ceil(rand(1)*nIndex);
    %choose random height
    randHeight = rand(1);

    %TO DO: Replace the ~= below with some other relationship:
    % geater than, less than, or equal?
    %if height is ?? compared to the probability value at this point in the
    %histogram, then select
    if (randHeight < probDist(r))
        break;
    end;    
end;
