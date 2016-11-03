function r=practicalLogReg2

%The goal of this part of the practical is to implement steepest descent
%and Newton's method on a 2 dimensional function

%This code finds the optimum of the Rosenbrocks function - the definition is in
%the routine Rosenbrock - look it up on the web!

%To make the routine work, you should calculate the first and second derivatives of this
%function by hand and fill in the routines below

%Things to investigate once you have gotten this to work
%1.  Change the routine to use steepest descent (Eq B.4 in the book's Appendix) rather than Newton's
%method
%2.  Implement the first derivative using finite differences rather than
%explicitly writing down the derivative (replace your implementaiton of RosenbrockDeriv)



%close all figures
close all;

%compute function on a grid to draw it

xIndex = -1.5:0.01:1.5;
wIndex = -1.5:0.01:1.5;
[x1 x2] = meshgrid(xIndex,wIndex);
x = [x1(:) x2(:)]';
w = Rosenbrock(x);
w = reshape(w,length(wIndex),length(xIndex));

%visualize range
figure; set(gcf,'Color',[1 1 1]);
imagesc(log(w+1));colormap(hot);
xlabel('x_{1}');ylabel('x_{2}');
set(gca,'YDir','Normal');
set(gca,'XTick',[2 length(xIndex)]);
set(gca,'YTick',[2 length(wIndex)]);
set(gca,'XTickLabel',{'-1.5','1.5'});
set(gca,'YTickLabel',{'-1.5','1.5'});

%TO DO complete this routine (see TODO's of the arguments passed to
%optimize(). 
startPosn = [-1.2;1];
tol = 10e-5;     %tolerance for finishing
[xPositions w] =optimize(startPosn,tol,@Rosenbrock,@RosenbrockDeriv,@Rosenbrock2ndDeriv);

%plot final result
hold on;
plot((xPositions(1,:)-xIndex(1))/(xIndex(2)-xIndex(1)),(xPositions(2,:)-wIndex(2))/(wIndex(2)-wIndex(1)),'g-');
plot((xPositions(1,end)-xIndex(1))/(xIndex(2)-xIndex(1)),(xPositions(2,end)-wIndex(1))/(wIndex(2)-wIndex(1)),'bo');

%plot actual minimum
plot((1-xIndex(1))/(xIndex(2)-xIndex(1)),(1-wIndex(1))/(wIndex(2)-wIndex(1)),'m+');

%==========================================================================

%Main optimization routine
function [xPositions w] =optimize(startPosn,tol,myFunction, myFunctionDeriv, myFunction2ndDeriv)

xPositions = startPosn;
w = myFunction(startPosn);

cIter = 1;
while(1)
    %compute derivative and Hessian at current positions
    derivVector = myFunctionDeriv(startPosn);
    hessian = myFunction2ndDeriv(startPosn);
    %compute update direction
    updateDirection = -inv(hessian)*derivVector;
    %do line search
    endPosn = lineSearch(startPosn,updateDirection,myFunction,tol);
    %store the end position and the best value so far
    xPositions = [xPositions endPosn];
    w = [w myFunction(endPosn)];
    %if didn't move very much then probably at maximum
    if(all(abs(endPosn-startPosn)<tol))
        break;  % Tip: may want to look up the all(.) function.
    else
        %new start position
        startPosn = endPosn;
    end;
    fprintf('Iter %d, Fn = %3.3f\n',cIter, myFunction(endPosn));
    cIter = cIter+1;
end;

%==========================================================================

function minX = lineSearch(startPosn,updateDirection,myFunction,tol)

%set options
options = optimset('Display','Off','TolX',tol);
startOffset = 0.0;

%run optimization
offset = fminsearch(@(offset)lineSearchFunction(offset,myFunction,updateDirection,startPosn),startOffset);
minX = startPosn+updateDirection*offset;


%==========================================================================

function r = lineSearchFunction(offset,myFunction,updateDirection,startPosn)

r = myFunction(startPosn+offset*updateDirection);



%==========================================================================
%Create function to optimize
function r= Rosenbrock(x)

%This is the Rosenbrock function (look it up!)
%It is a simple polynomial equation but it is quite hard to find the exact
%minimum!
x1 = x(1,:);
x2 = x(2,:);
r = 100*(x2-x1.^2).^2+(1-x1).^2;



%==========================================================================
function deriv= RosenbrockDeriv(x)
%return partial derivatives of Rosenbrock function (a 2x1 vector)
% and note that the input argument, x, is 2D.

%TO DO  - compute the first derivatives of Rosenbrock's function
%Replace this:
deriv = ones(2,1);

if ~(size(deriv, 1) == 2 && size(deriv, 2) == 1)
    error('deriv should be a 2 x 1 vector');
end


%==========================================================================
%Compute 2nd derivative of function
function hessian= Rosenbrock2ndDeriv(x)
%return hessian of Rosenbrock function (a 2x2 matrix)

%TO DO  - compute the second derivatives of Rosenbrock's function
%Replace this:
hessian = eye(2,2);

if ~(size(hessian, 1) == 2 && size(hessian, 2) == 2)
    error('hessian should be a 2 x 2 matrix');
end
