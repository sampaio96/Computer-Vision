function practicalRegress4
%Goal: to perform Gaussian process regression of pose data.

% The implementation of Gaussian process regression 
% regresses images of faces against poses. The code takes about a minute 
% to run as the dimensionality and number of training examples are both quite high.

%Look through the code and make sure you understand everything it is doing
%Head pose estimation is a difficult problem - no-one really knows how to solve this yet
%This method is not quite cutting edge - we can get the
%mean average error down to ~10 degrees with a different approach.
%Even humans cannot do this very well!

%close all previous figures
close all;

%load in training data
load('PoseRegressionData.mat','dataTrain','dataTest','dataTrainPP','dataTestPP','poseTrain','poseTest');

%find size of preprocssed image
[imY imX I] = size(dataTrainPP);
%reshape training data into columns of matrix
dataTrainPP = reshape(dataTrainPP,[],I);

%place one at the top of every data example
%dataTrainPP = [ones(1,I);dataTrainPP];
[D I] = size(dataTrainPP);

%define hyperparameter
sigmaPSq = 1000000;

%The goal of this routine is to take data x and w and the prior variance
sigmaSq = estNoiseParameter(dataTrainPP,poseTrain,sigmaPSq);

%now fit model
W = inv(kernel(dataTrainPP,dataTrainPP)+(sigmaSq/sigmaPSq)*eye(I));
WKXXw = W*kernel(dataTrainPP,dataTrainPP)*poseTrain;
%find size of test data
[imY imX nTestData] = size(dataTestPP);
%reshape training data into columns of matrix
dataTestPP = reshape(dataTestPP,[],nTestData);

%for each test data, get mean and variance of prediction
poseTestPredictMean = zeros(nTestData,1);
poseTestPredictVar = zeros(nTestData,1);
for(cTestData = 1:nTestData)
    thisData = dataTestPP(:,cTestData);
    poseTestPredictMean(cTestData) = (sigmaPSq/sigmaSq)*kernel(thisData,dataTrainPP)*poseTrain-(sigmaPSq/sigmaSq)*kernel(thisData,dataTrainPP)*WKXXw;
    poseTestPredictVar(cTestData) = sigmaPSq*kernel(thisData,thisData)-sigmaPSq*kernel(thisData,dataTrainPP)*W*kernel(dataTrainPP,thisData)+sigmaSq;
end;

%draw figure of predictions vs. ground truth
figure; set(gcf,'Color', [1 1 1]);
plot(poseTest,poseTestPredictMean,'r.');
xlim([-90 90]); xlabel('Actual Pose');
ylim([-90 90]); ylabel('Predicted Pose');
hold on; 
plot([-90 90],[-90 90],'k-');
axis square;

%compute statistics of how well we have done
covMat = cov([poseTest poseTestPredictMean]);
PPMCC = covMat(2)/sqrt(covMat(1)*covMat(4));
MAE = mean(abs(poseTest-poseTestPredictMean));
fprintf('Pearson product moment coefficient = %f\n',PPMCC);
fprintf('Mean average error = %f\n',MAE);

%draw some of the data
figure; set(gcf,'Color',[1 1 1]);
set(gcf,'Position',[    47    37   898   654]);
nTestData = size(dataTest,4);
randOrder = randperm(nTestData);
for (cData = 1:16)
    subplot(4,4,cData);
    imshow(uint8(dataTest(:,:,:,randOrder(cData))));
    title(sprintf('Est: %2.1f, True: %2.1f',poseTestPredictMean(randOrder(cData)),poseTest(randOrder(cData))));
end;



%==========================================================================
%==========================================================================

%returns negative log marginal likelihood of data
%i.e. marginal likelihood is likelihood after marginalizing over phi
function L = negLogMarginalLike(sigmaSq,kernelXX,w,sigmaPSq)

%take exponential as still in log form
sigmaSq = exp(sigmaSq);

%compute mean and variance parameters
muParam = zeros(length(w),1);
covarParam = kernelXX*sigmaPSq+sigmaSq*eye(size(kernelXX,2));
%return negative log marginal
L = -getLogGaussianLike(w,muParam,covarParam);


%==========================================================================
%==========================================================================


%returns log of normal pdf
function L = getLogGaussianLike(w,muParam,covarParam)

D = length(w);
L =-0.5*D*log(2*pi)-0.5 *logDet(covarParam)-0.5*(muParam-w)'*inv(covarParam)*(muParam-w);


%==========================================================================
%==========================================================================


%returns log of determinant of matrix efficiently
%The determinant is a very small number so we cannot compute it and then
%take the log as matlab can't represent it.  So this is computed in a
%sneaky way!
function ld = logDet(A)

[U L V] = svd(A);
ld = sum(log(diag(L)));



%==========================================================================
%==========================================================================

%returns RBF kernel matrix
function K = kernel(X1,X2);

lengthScale=3000;

I1 = size(X1,2);
I2 = size(X2,2);

%create each row of kernel matrix separately.
K = zeros(I1,I2);
for (c1 = 1:I1)
   %compute distance between this example and all other examples
   diff = sum((repmat(X1(:,c1),1,I2)-X2).^2)/lengthScale ;   
   %store in kernel matrix
   K(c1,:) = exp(-diff);
end

%==========================================================================
%==========================================================================


function sigmaSq = estNoiseParameter(dataTrainPP,poseTrain,sigmaPSq);
%initial estimate for variance
sigmaSqInit = std(poseTrain);
%take logarithm (so optimizer doesn`t have to constrain to be positive)
sigmaSqInit = log(sigmaSqInit);
%precompute kernel so that it doesn't have to repeatedly do this
kernelXX = kernel(dataTrainPP,dataTrainPP);
%fit variance of data - this routine finds the value of sigmaSq that
%minimizes the negative log marginal likelihood
sigmaSq = fminsearch(@(sigmaSq) negLogMarginalLike(sigmaSq, kernelXX,poseTrain,sigmaPSq),sigmaSqInit);
%exponentiate to return to positive value
sigmaSq = exp(sigmaSq)