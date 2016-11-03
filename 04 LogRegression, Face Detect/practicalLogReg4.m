function r=practicalLogReg4
%The goal of this part of the practical is to implement logistic
%regression for classifying faces and non-faces

%This is a complete working bit of code. Your goals are to 
%(i) Look at the code and understand it
%(ii) Investigate what happens as you increase the amount of training data
%- try 750, 1000, 1500,2000,2000,3000, 4000 examples.  Does it generalize better?
%(iii) Try learning with gradient descent with 4000 examples. What happens?
%(iii) Convert this to a non-linear logistic regression algorithm by
%transforming the data before running the routine (don't forget to
%transform the test data as well). Write a routine to transform each
%data point to a 500x1  vector by evaluating it against 500 radial basis
%functions.  The centers of these functions can be the first 500 data
%points.  You should experiment with the standard deviation, but somewhere
%in the range 1-100 should be a good start.

%close all figures
close all;
%set seeds so always get same random numbers
randn('seed',2);
rand('seed',2);

%load in training data
load('FaceDetectData.mat','x','y','xTest','yTest');
%select the amount of training data
nTrainData = 600;
x = x(:,:,1:nTrainData);
y = y(1:nTrainData);
%draw some of the data examples
figure; set(gcf,'Color',[1 1 1]);
for (cFig = 1:16);
    subplot(4,4,cFig);
    imagesc(x(:,:,cFig)); axis off; axis image;colormap(gray);
    if(y(cFig))
        title('Face');
    else
        title('Non-Face');
    end;
end;

%vectorize data
xVec = reshape(x,size(x,1)*size(x,2),size(x,3));
%prepend ones
xVec = [ones(1,size(xVec,2));xVec];
%do optimization
startPosn = [1; 0.0001*randn(size(xVec,1)-1,1)];
tol = 10e-8;
[phiValues yPredicted] =optimize(startPosn,tol,@NegLogProb,@NegLogProbDeriv,@NegLogProbHessian,xVec,y);

phi = phiValues(:,end);

%get summary statistics
prob = sig(phi'*xVec);
class = prob>0.5;
nCorrect = sum(class==y);
fprintf('Training Data: Classified %3.3f percent correct\n',100*nCorrect/length(y));

%compute for test data
xTestVec = reshape(xTest,size(xTest,1)*size(xTest,2),size(xTest,3));
xTestVec = [ones(1,size(xTestVec,2));xTestVec];
%get summary statistics
prob = sig(phi'*xTestVec);
class = prob>0.5;
nCorrect = sum(class==yTest);
fprintf('Test Data: Classified %3.3f percent correct\n',100*nCorrect/length(yTest));



%draw some of the test examples
figure; set(gcf,'Color',[1 1 1]);
for (cFig = 1:16);
    subplot(4,4,cFig);
    imagesc(xTest(:,:,cFig)); axis off; axis image;colormap(gray);
    if(prob(cFig)>0.5)
        title('Face');
    else
        title('Non-Face');
    end;
end;


%==========================================================================

%Main optimization routine
function [xPositions functionValues] =optimize(startPosn,tol,myFunction, myFunctionDeriv, myFunction2ndDeriv,x,y)

xPositions = startPosn;
functionValues = myFunction(startPosn,x,y);

fprintf('Iter %d, Fn = %3.3f\n',0, myFunction(startPosn,x,y));

cIter = 1;
MAX_ITER = 50;
while(cIter<MAX_ITER)
    %compute derivative and Hessian at current positions
    derivVector = myFunctionDeriv(startPosn,x,y);
    hessian = myFunction2ndDeriv(startPosn,x,y);
    %compute update direction (based on gradient descent) 
    updateDirection = -hessian\derivVector;
    %do line search
    endPosn = lineSearch(startPosn,updateDirection,myFunction,tol,x,y);
    %store the end position and the best value so far
    xPositions = [xPositions endPosn];
    functionValues = [functionValues myFunction(endPosn,x,y)];
    %Print out statistics
    prob = sig(endPosn'*x);
    class = prob>0.5;
    nCorrect = sum(class==y);
    fprintf('Iter %d, Fn = %3.3f, percentCorrect = %3.3f\n',cIter, myFunction(endPosn,x,y),100*nCorrect/length(y));

    %if didn't move very much then probably at maximum
    if(min(abs(endPosn-startPosn)<tol)|functionValues(end)<10.0)
        break;
    else
        %new start position
        startPosn = endPosn;
    end;
    cIter = cIter+1;
end;

%==========================================================================

function minX = lineSearch(startPosn,updateDirection,myFunction,tol,x,y)

%set options
options = optimset('Display','off','TolX',tol,'TolFun',tol);
startOffset = 0.0;

%run optimization
offset = fminsearch(@(offset)lineSearchFunction(offset,myFunction,updateDirection,startPosn,x,y),startOffset,options);
minX = startPosn+updateDirection*offset;


%==========================================================================

function r = lineSearchFunction(offset,myFunction,updateDirection,startPosn,x,y)

r = myFunction(startPosn+offset*updateDirection,x,y);

%==========================================================================
function r = sig(a)

r = 1./(1+exp(-a));
r = r;


%==========================================================================
%returns log probability of training data (a scalar)
function L= NegLogProb(phi, x,y)

smallNum = 10e-15;
%added a small number to stop getting log(0);
L = y.*log(sig(phi'*x)+smallNum)+(1-y).*log(1-sig(phi'*x)+smallNum);
L = -sum(L);
L = L;

%==========================================================================
%Compute derivative of log probability (3x1 vector)
function deriv= NegLogProbDeriv(phi,x,y)
%return derivative of Rosenbrock function (a 2x1 vector)

deriv = ones(3,1);
deriv  = x.*repmat(sig(phi'*x)-y,size(x,1),1);
deriv = -sum(deriv,2);


%==========================================================================
%Compute 2nd derivative of log probability (3x3 matrix)
function hessian= NegLogProbHessian(phi,x,y)

hessian = ones(3,3);
weights = sig(phi'*x).*(1-sig(phi'*x));
weightedX = x.*repmat(weights,size(x,1),1);
hessian = -weightedX*x';
