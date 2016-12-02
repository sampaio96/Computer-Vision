function r =  practical6b
%In practical 6b, we experiment with sampling from
%probability distributions using MCMC methods. 


%close all previous figures
close all;

%First we define a the mean and covariance of a 2D Gaussian
gaussMean = [1.0000; 2.0000];
gaussCov =  [1.6370 0.4210;...
             0.4210 0.1598];
         
%draw the Gaussian
figure; set(gcf,'Color',[ 1 1 1]);
drawGaussianOutline(gaussMean,gaussCov); hold on;
plot(gaussMean(1),gaussMean(2),'r.');
set(gca,'Box','Off'); 
xlim([-5 5]); ylim([-5 5]);

%initialize first point in chain
position = [0;0];
plot(position(1),position(2),'b.','MarkerSize',5);

%define number of samples
nSample = 1000;

%loop over number of samples
for (cSample = 1:nSample)
        %TO DO:  calculate mean and variance of conditional distribution in
        %x direction given current y positition - look at the notes for the
        %Gaussian distribution. Replace this:
        gaussMeanXGivenY = gaussMean(1) + gaussCov(1,2)*inv(gaussCov(2,2))*(position(2)-gaussMean(2));
        gaussVarXGivenY = gaussCov(2,2)-gaussCov(1,2)'*inv(gaussCov(1,1))*gaussCov(1,2);
        
        %TO DO:  draw a sample from this distribution and update x!!!! - use
        %randn to do this. Replace this:
        newPosition(1) = gaussMeanXGivenY+gaussVarXGivenY*randn(1);
        
        %copy the y position
        newPosition(2) = position(2);
        
        %display this update
%         plot([position(1) newPosition(1)],[position(2) newPosition(2)],'b-');
        
        %update the original position
        position = newPosition;
        
        %TO DO:  calculate mean and variance of conditional distribution in
        %y direction given current x positition - look at the notes for the
        %Gaussian distribution. Replace this:
        gaussMeanYGivenX = gaussMean(2) + gaussCov(2,1)*inv(gaussCov(1,1))*(position(1)-gaussMean(1));
        gaussVarYGivenX = gaussCov(1,1)-gaussCov(2,1)'*inv(gaussCov(2,2))*gaussCov(2,1);
        
        %TO DO:  draw a sample from this distribution and update y - use
        %randn to do this. Replace this
        newPosition(2) = gaussMeanYGivenX+gaussVarYGivenX*randn(1);
        
        %copy the y position
        newPosition(1) = position(1);

        %display this update
%         plot([position(1) newPosition(1)],[position(2) newPosition(2)],'b-');
        
        %update the original position
        position = newPosition;
        
        %display the new point
        plot(position(1),position(2),'b.');
end;

%TO DO now remove the two lines of code that plot the path between points and use
%this method to generate 1000 samples from the Gaussian

%*************************************************************************
%Please note that this is not the usual way to sample from an n-d Gaussian!  
%A much more efficient way to do this is to take:

%newSample = gaussMean+chol(gaussCov)*randn(2,1);

%the routine "chol" is the cholesky decomposition - a matrix square root
% so that A = chol(A)'*chol(A);
%**************************************************************************

%=================================================================== 
%===================================================================

%draw 2DGaussian
function r= drawGaussianOutline(m,s)

hold on;
angleInc = 0.01;

c = [1.0 0.0 0.0];

for (cAngle = 0:angleInc:2*pi)
    angle1 = cAngle;
    angle2 = cAngle+angleInc;
    [x1 y1] = getGaussian2SD(m,s,angle1);
    [x2 y2] = getGaussian2SD(m,s,angle2);
    plot([x1 x2],[y1 y2],'k-','LineWidth',2,'Color',c);
end

%===================================================================
%===================================================================

%find position of in xy co-ordinates at 2SD out for a certain angle
function [x,y]= getGaussian2SD(m,s,angle1)

if (size(s,2)==1)
    s = diag(s);
end;

vec = [cos(angle1) sin(angle1)];
factor = 4/(vec*inv(s)*vec');

x = cos(angle1) *sqrt(factor);
y = sin(angle1) *sqrt(factor);

x = x+m(1);
y = y+m(2);

