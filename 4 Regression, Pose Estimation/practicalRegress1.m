function r=practicalRegress1

%The goal of this practical is to investigate methods for regression
%In part 1 we look at linear regression with maximum likelihood learning
%In part 2 we look at linear regression with Bayesian learning
%In part 3 we look at non-linear regression
%In part 4, observe Gaussian process regression; step through in debugger.


%Work your way through these Matlab examples filling in code where it says
%"TO DO"


%close all figures
close all;

%set seeds for random number generators 
%so we get the same random numbers each time
rand('seed',5);
randn('seed',5);


%define actual parameters
offsetActual = 1.5;  % this is phi_0 in the notes
slopeActual = -0.8;  % this is phi_1 in the notes
noiseActual = 0.01;  % this is sigma^2 in the notes


%generate some data 
nData = 10;
x = rand(1,nData);
w = offsetActual + slopeActual * x + sqrt(noiseActual) * randn(1,nData);
%Transpose w so that it is a column vector like in the notes
w = w';


%display this data
figure; set(gcf,'Color',[1 1 1]);
plot(x,w','r.');
xlabel('Data, x'); ylabel('World, w');
set(gca,'Box','Off'); 
xlim([0 1]);ylim([0 2]);

%now the problem is to fit the data
%TODO  - fill in this routine (at the bottom of this file, where its
%skeleton can already be found)
[offsetEst slopeEst noiseEst] = fitLinearRegressionML(x,w);

%Let's draw a figure to show the posterior probability 
%of the world as a function of the data
%For every x, the model predicts a normal distribution over w

x2 = 0:0.005:1;
% TO DO Inference: Compute the mean and variance of the prediction for each element of
% x2.  Replace this:
x2a = [ones(size(x2)); x2];
predMean = (x2a.' * [offsetEst; slopeEst]).';
predVar = noiseEst * ones(size(x2));


% performing checks on the size of the data created 
if ~(isvector(predMean) && size(predMean, 2)==length(x2))
    error('predMean should be same size as variable x2');
elseif ~(isvector(predVar) && size(predVar, 2)==length(x2))
    error('predVar should be same size as variable x2');
end

%compute the pixel colours for the image - one Gaussian in each column
wFig = 0:0.005:2;
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


%The goal of this routine is to take data x and w and fit the three
%parameters of the linear regression model
function [offsetEst slopeEst noiseEst] = fitLinearRegressionML(x,w)

%number of data points 
nData = size(x,2); % I

%TODO add a one to the start of each data example x

x = [ ones(1,nData); x ];

%TODO compute phi vector (contains offset and slope)
%REPLACE THIS:
phi = ones(1,2);
phi = inv(x * x.') * x * w;

%extract the slope and offset from this vector
offsetEst = phi(1);
slopeEst = phi(2);

%compute the variance parameter
%Replace this
noiseEst = (w - x.' * phi).' * (w - x.' * phi) / nData;

% performing checks on the data created
if ~(isvector(phi) && length(phi)==2)
    error('phi should have size 2 x 1');
elseif ~isscalar(noiseEst)
    error('noiseEst should be scalar');
end
