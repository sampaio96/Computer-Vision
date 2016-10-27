function r=practicalRegress3

%The goal of this practical is to implement non-linear regression method

%Work your way through these Matlab examples filling in code where it says
%"TO DO"

%TO DO 
%First Run this file.  It implements linear regression on a dataset that is
%clearly not very linear.  It's a horrible fit!
%Your job is to improve this by doing non-linear regression (as in Figure
%8.6).


%close all figures
close all;

load('nonLinearData.mat','x','w');
%keep a copy of the original data
[nDim nData] = size(x);

%display this data
figure; set(gcf,'Color',[1 1 1]);
plot(x,w','r.');
xlabel('Data, x'); ylabel('World, w');
set(gca,'Box','Off'); 
xlim([-3 3]);ylim([-1.5 1.5]);

%TO DO 
%To implement non-linear regression, we simply pass the data through
%a non-linearity before we do the regression
%Fill in the routine (code below)
xTransform = nonLinearTransform(x);

%now the problem is to fit the data
[offsetEst slopeEst noiseEst] = fitLinearRegressionML(xTransform,w);

%Let's make some draw a figure to show the posterior probability 
%of the world as a function of the data
%For every x, the model predicts a normal distribution over w

x2 = -3:0.01:3;
x2Transform = nonLinearTransform(x2);
% Compute the mean and variance of the prediction for each element of
% x2.  
predMean =slopeEst'*x2Transform+offsetEst;
predVar =noiseEst*ones(size(x2Transform,2));

%compute the pixel colours for the image - one Gaussian in each column
wFig = -1.5:0.01:1.5;
postFun = zeros(length(wFig),length(x2));
for (cX= 1:length(x2))
    postFun(:,cX)= (1./sqrt(2*pi*predVar(cX))).*exp(-0.5 *((wFig-predMean(cX)).^2)./predVar(cX));
end;

%draw the figure
figure; set(gcf,'Color',[1 1 1]);
imagesc(postFun); colormap(hot);
set(gca,'YDir','Normal');
hold on;
%draw the points on top - they need to be rescaled so that they are in the 
%units of the pixels used to draw the picture
plot((x-x2(1))/(x2(2)-x2(1)), (w'-wFig(1))/(wFig(2)-wFig(1)),'bo');
%draw the units on the graph
set(gca,'XTick',[1 length(x2)]);
set(gca,'XTickLabel',{'0','1'});
set(gca,'YTick',[1 length(wFig)]);
set(gca,'YTickLabel',{'0','1'});
xlabel('x');
ylabel('w');
title('Pr(w|x)');

%TO DO
%When you've got this working you could try change this to work with the 
%arc tan functions if you like (figure 8.7)
% Note: Try some parameters, see if you can find better ones by trial and
% error. Wouldn't it be nice to have a huge (or infinite?) set of basis
% functions? Possible with Gaussian Process Regression!




%The goal of this routine is to take data x and w and fit the three
%parameters of the linear regression model
function [offsetEst slopeEst noiseEst] = fitLinearRegressionML(x,w)

%number of data points 
nData = size(x,2);

%add a one to the start of each data example x
x = [ones(1,nData);x]; % x : 1+nDim x nData

%compute phi vector (contains offset and slope)
phi = inv(x*x')*x*w;

%extract the slope and offset from this vector
offsetEst = phi(1);
slopeEst = phi(2:end);

%compute the variance parameter
noiseEst = (w-x'*phi)'*(w-x'*phi)/nData;


%The goal of this method is to pass the data through a non-linear
%transformation
function x2= nonLinearTransform(x);

%TO DO 
%FIRST REMOVE THESE TWO LINES
% x2 = x;
% return;


nData= size(x,2);
%create space for output data
x2 = zeros(6,nData);
%For each data points
for (cData = 1:nData)
    thisX = x(:,cData);
    %TO DO: perform non-linear transformation on thisX:
    %To do this evaluate this X against the 6 RBF functions in Figure 8.6b 
    %They are 6 Gaussians with means
    % -2.5 -1.5 -0.5 0.5 1.5 and 2.5, and variances of 0.16;
    %replace this:
    var = 0.16;
    means = [-2.5 -1.5 -0.5 0.5 1.5 2.5];
    transformedX = zeros(6,1);
    transformedX = exp(-(thisX-(means)).^2/(2*var));
    
    %store
    x2(:,cData) = transformedX;
end;

% checking for size of x2
if ~(all(size(x2) == [6, nData]))
    error('x2 should have dimensions 6 x nData')
end