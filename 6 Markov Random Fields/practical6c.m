function r=practical6c
%The goal of this practical is to sample from a binary Markov random field
%Once again we are going to use Gibbs sampling.  We select one dimension
%(pixel) at a time, and consider the conditional distribution of this pixel
%given all of the others.  Then we sample from this distribution and
%repeat.

%close all prevoius figures
close all;

%define costs of Markov random field -costs for neighbouring labels.
MRFCosts = [0.001 0.60;...   
            0.60 0.001];

        
        
%define size of label field
imX = 20; imY = 20;

%choose initial random labels for field
label = round(rand(imY,imX));

%display initial random labels
figure; set(gcf,'Color',[1 1 1]);
imagesc(label); axis off; axis image; colormap(gray);

%define number of iterations
nIter = 1000;

%run through iterations
for (cIter = 1:nIter)
    %define random order of pixels to update
    pixelOrder = randperm(imX*imY);
    [pixelOrderX pixelOrderY]=ind2sub([imY imX],pixelOrder);
    %run through each pixel
    for (cPixel = 1:imY*imX)
        %choose which pixel to update
        thisX= pixelOrderX(cPixel);
        thisY= pixelOrderY(cPixel);
        %set this label to zero
        label(thisY,thisX) = 0;
        %calculate the probability of the MRF 
        %TO DO - fill in this routine (below)
        prLabelEquals0 = calcMRFProbability(label,MRFCosts);
        %set this label to one
        label(thisY,thisX) = 1;
        %calculate the probability of the MRF
        prLabelEquals1 = calcMRFProbability(label,MRFCosts);
        %now we will sample from conditional probability distribution        
        %TO DO normalize the two probabilities so that they sum to one
        Z = prLabelEquals0+prLabelEquals1;
        prLabelEquals0 = prLabelEquals0/Z;
        prLabelEquals1 = prLabelEquals1/Z;
        %TO DO  Generate a random number between 0 and 1. If it is
        %less than prLabelEquals0, then set the label to zero, otherwise set
        %it to one. Replace this:
        label(thisY,thisX) = double(rand(1)<prLabelEquals0);
    
        %draw image 
        imagesc(label); axis off; axis image; drawnow;   
    end;
end;

%if you have got this right, the image should gradually become smooth.  Over time, 
%the images are samples from the imY*imX probability distribution over
%label fields.

%TO DO:  Now look at the computation - almost all of the terms that you compute
%are the same for both prLabelEquals1 and prLabelequals0 and hence they
%cancel when you normalize.  How would you re-implement this so that you
%did not do all of this unnecessary compution?

    
%==========================================================================
%==========================================================================

%the goal of this routine is to calculate the probability of 
%the Markov random field (up to the unknown scaling factor!)
function prMRF = calcMRFProbability(label,MRFCosts)

%find size of image
[imY imX] = size(label);

%initialize total cost - sum of psi terms
U = 0;
Z = zeros(imY, imX);

%for each pixel
for(cPixelY = 1:imY)
    for (cPixelX = 1:imX)
        %add cost for neighbour above this pixel
        %TO DO - fill in routine get cost
        if(cPixelY>1)
            U = U + getCost(label(cPixelY,cPixelX),label(cPixelY-1,cPixelX),MRFCosts);
            Z(cPixelY, cPixelX) = Z(cPixelY, cPixelX) + getCost(label(cPixelY,cPixelX),label(cPixelY-1,cPixelX),MRFCosts);
        end;
        %add cost for neighbour below this pixel
        if(cPixelY<imY)
            U = U + getCost(label(cPixelY,cPixelX),label(cPixelY+1,cPixelX),MRFCosts);
            Z(cPixelY, cPixelX) = Z(cPixelY, cPixelX) + getCost(label(cPixelY,cPixelX),label(cPixelY+1,cPixelX),MRFCosts);
        end;
        %add cost for neighbour to the left
        if(cPixelX>1)
            U = U + getCost(label(cPixelY,cPixelX),label(cPixelY,cPixelX-1),MRFCosts);
            Z(cPixelY, cPixelX) = Z(cPixelY, cPixelX) + getCost(label(cPixelY,cPixelX),label(cPixelY,cPixelX-1),MRFCosts);
        end;
        %add cost for neighbour to the right(!!!) this pixel
        if(cPixelX<imX)
            U = U + getCost(label(cPixelY,cPixelX),label(cPixelY,cPixelX+1),MRFCosts);
            Z(cPixelY, cPixelX) = Z(cPixelY, cPixelX) + getCost(label(cPixelY,cPixelX),label(cPixelY,cPixelX+1),MRFCosts);
        end;
    end;
end;

%TO DO - calculate cost of MRF. Replace this:
Z = sum(sum(Z));
prMRF = exp(-U)/Z;

%==========================================================================
%==========================================================================

%The goal of this routine is to return the cost of the MRF
function cost = getCost(label1,label2,MRFCosts)
%TO DO - fill in this routine - Replace this
cost = -log(MRFCosts(label1+1, label2+1));

    
    
    
    
    
    