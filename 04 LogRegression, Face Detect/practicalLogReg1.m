function r=practicalLogReg1

%The goal of this practical is to investigate the logistic regression model
%In part 1 we investigate line search 
%In part 2 we will build models to optimize a function
%In part 3 we will fit the logistic regression model to some data
%In part 4 we will apply this to the face detection problem (inspect this!)
% Fill in sections labeled "TO DO".

%close all figures
close all;

%Compute function over some range. Can think of myFunction as a "black box".
x = 0:0.01:1;
w = zeros(size(x));
for (cX = 1:length(x))
    w(cX) = myFunction(x(cX));
end;

%visualize range
figure; set(gcf,'Color',[1 1 1]);
plot(x,w,'r-');
xlabel('x');ylabel('f(x)');
ylim([-1.5 0.5]);

%define parameters for line search
startSearch = 0; %beginning of search region
endSearch = 1;   %end of search region
tol = 10e-5;     %tolerance for finishing
%TO DO complete this routine
[minX minValue] = lineSearch(startSearch,endSearch,tol,@myFunction);
% Note: '@myFunction' means we are passing in a function-handle. This way,
% lineSearch() is receiving the function myFunction() as one of the
% parameters passed to it.


%plot final result
hold on;
plot(minX,minValue,'bo');

%Note that there is a built in MATLAB routine that will do this for you
%set options
options = optimset('Display','On','TolX',tol);
startPosition = 0.0;

%run optimization. TO DO: Look up HELP for this function.
minX = fminsearch(@myFunction,startPosition);

%draw result
minValue = myFunction(minX);
plot(minX,minValue,'g+');

%
% Now repeat the process using a different "black box": myFunction2!
%
w = zeros(size(x));
for (cX = 1:length(x))
    w(cX) = myFunction2(x(cX));
end;
% TO DO: run your linesearch on this new training data, trying
% different initializations.


%==========================================================================

%main line search function
function [minX minValue] = lineSearch(startSearch,endSearch,tol,optFunction)

%define points a,b,c,d
a= startSearch;
d = endSearch;
while((d-a)>tol)
    fprintf('Current Search Region is %6.6f to %6.6f\n',a,d);
    %TO DO:  calculate intermediate positions b and c, 1/3 and 2/3 of the way 
    % through interval (between a and d), respectively
    
    %replace this
    break;
    
    
    %TO DO evaluate function at points b and c
    
    
    %TO DO:  update point a or d depending on the values of those evaluations
    
end;

%calculate value at minimum
minX = a;
minValue = optFunction(a);




%==========================================================================
%Create function to optimize
function r= myFunction(x)

if (x<0)
    r = -x;
elseif (x>1)
   r =  -sin((1*1.8).^2)+x-1;
else
    r = -sin((x*1.8).^2);
end;

%==========================================================================
%Create function to optimize
function r= myFunction2(x)

if (x<0)
    r = -x-exp(-0.5*((-0.2)/0.025).^2);;
elseif (x>1)
   r = -sin((1*1.8).^2)+x-1;
else
    r = -sin((x*1.8).^2);
    r = r-exp(-0.5*((x-0.2)/0.025).^2);
end;