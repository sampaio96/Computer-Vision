function r=practicalLogReg3

%The goal of this part of the practical is to implement Newton's method
%for fitting logistic regression with 2D data

%As for the last part, we have implemented the main function, so you should
%implement the first and second derivatives.

%When you have this working, compare the number of iterations for steepest
%descent vs. Newton's method

%close all figures
close all;

%set seed so always get same random numbers
rand('seed',3)
randn('seed',3)

%number of data
nData = 20;

%generate some random values from Gaussian distribution 
x1 = randn(2,nData)*0.15+0.65; t1 = ones(1,nData);
x2 = randn(2,nData)*0.19+0.4; t2 = zeros(1,nData);

%show data before
figure; set(gcf,'Color',[1 1 1]);
plot(x1(1,:),x1(2,:),'go');hold on
plot(x2(1,:),x2(2,:),'bo');
axis off; axis square;

%concatenate data
x = [ones(1,length(x1)+length(x2));x1 x2];
w = [t1 t2];



startPosn = [1;-2.0;2.0];
tol = 10e-5;
[phiValues wPredicted] =optimize(startPosn,tol,@NegLogProb,@NegLogProbDeriv,@NegLogProbHessian,x,w);


%compute function on a grid to draw it
xIndex = -0.2:0.01:1.2;
wIndex = -0.2:0.01:1.2;
[x1Mesh x2Mesh] = meshgrid(xIndex,wIndex);
x = [x1Mesh(:) x2Mesh(:)]';
x = [ones(1,size(x,2)); x];
phi = phiValues(:,end);
prob = sig(phi'*x);


prob = reshape(prob,length(wIndex),length(xIndex));

%visualize range
figure; set(gcf,'Color',[1 1 1]);
imagesc(prob);colormap(hot);
xlabel('x_{1}');ylabel('x_{2}');
set(gca,'YDir','Normal');
set(gca,'XTick',[2 length(xIndex)]);
set(gca,'YTick',[2 length(wIndex)]);
set(gca,'XTickLabel',{'-0.2','1.2'});
set(gca,'YTickLabel',{'-0.2','1.2'});

%plot the data on top
hold on;
plot((x1(1,:)-xIndex(1))/(xIndex(2)-xIndex(1)),(x1(2,:)-wIndex(2))/(wIndex(2)-wIndex(1)),'go');
plot((x2(1,:)-xIndex(1))/(xIndex(2)-xIndex(1)),(x2(2,:)-wIndex(2))/(wIndex(2)-wIndex(1)),'bo');


%==========================================================================

%Main optimization routine
function [xPositions functionValues] =optimize(startPosn,tol,myFunction, myFunctionDeriv, myFunction2ndDeriv,x,w)

xPositions = startPosn;
functionValues = myFunction(startPosn,x,w);

cIter = 1;
while(1)
    %compute derivative and Hessian at current positions
    derivVector = myFunctionDeriv(startPosn,x,w);
    hessian = myFunction2ndDeriv(startPosn,x,w);
    %compute update direction 
    updateDirection = -inv(hessian)*derivVector;
    %do line search
    endPosn = lineSearch(startPosn,updateDirection,myFunction,tol,x,w);
    %store the end position and the best value so far
    xPositions = [xPositions endPosn];
    functionValues = [functionValues myFunction(endPosn,x,w)];
    %if didn't move very much then probably at maximum
    if(min(abs(endPosn-startPosn)<tol))
        break;
    else
        %new start position
        startPosn = endPosn;
    end;
    fprintf('Iter %d, Fn = %3.3f\n',cIter, myFunction(endPosn,x,w));
    cIter = cIter+1;
end;

%==========================================================================

function minX = lineSearch(startPosn,updateDirection,myFunction,tol,x,w)

%set options
options = optimset('Display','Off','TolX',tol);
startOffset = 0.0;

%run optimization
offset = fminsearch(@(offset)lineSearchFunction(offset,myFunction,updateDirection,startPosn,x,w),startOffset);
minX = startPosn + updateDirection * offset;


%==========================================================================

function r = lineSearchFunction(offset,myFunction,updateDirection,startPosn,x,w)

r = myFunction(startPosn+offset*updateDirection,x,w);

%==========================================================================
function r = sig(a)

r = 1./(1+exp(-a));
r = r;


%==========================================================================
%Returns negative log probability of training data (so output a scalar)
function L= NegLogProb(phi, x,w)

smallNum = 10e-6;
%added a small number to stop getting log(0);
L = w.*log(sig(phi'*x)+smallNum)+(1-w).*log(1-sig(phi'*x)+smallNum);
L = -sum(L);

%==========================================================================
%Compute derivative of negative log probability (3x1 vector)
function deriv= NegLogProbDeriv(phi,x,w)
%TO DO - compute this derivative.
%replace this.
deriv = ones(3,1);

if ~(size(deriv, 1) == 3 && size(deriv, 2) == 1)
    error('deriv should be a 3 x 1 vector');
end

%==========================================================================
%Compute 2nd derivative of negative log probability (3x3 matrix)
function hessian= NegLogProbHessian(phi,x,w)

%TO DO - compute this matrix.
%Replace this:
hessian = eye(3,3);

if ~(size(hessian, 1) == 3 && size(hessian, 2) == 3)
    error('hessian should be a 3 x 3 matrix');
end