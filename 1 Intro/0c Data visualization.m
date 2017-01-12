%Practical 0c. Figures and display

%The goal of this practical is to investigate the plotting and display
%functions of Matlab.

%Close all previous figures and clear all memory
close all;
clear all;

%let's create some data
x = 1:20;
y = exp(x/10);
z = sin(x);

%we can create a new figure using
h = figure
% the command figure returns a handle h that will help us modify the figure

%Each figure has a set of properties associated with it

%TO DO
%retrieve the figures properties using the command get(h);
get(h)
%TO DO 
%retrieve the current colour of the figure by using the syntax
%get(h,'PropertyName')
get(h,'Color')
%TO DO 
%change the colour of the figure to white using the syntax
%set(h,'PropertyName',newValue)
set(h,'Color',[1 1 1])
%TO DO 
%create a new figure and store the handle in variable g
%change the position of and size of this figure using the 'Position'
%property
g = figure
set(g,'Position',[440 478 460 400])
%TO DO 
%select the first figure again by using
%figure(h)
figure(h)
%Let us plot a graph of x against y using the command
j = plot(x,y,'r>-');
%The 'r>.' is a line specification.  It
%means draw a red line that is continuous with triangular markers

%TO DO 
%Change this command to draw a dashed green line with circular markers. You
%can find out how to do this by looking up 'linespec' in the help.  Use the
%command 'doc linespec' to do this.
j = plot(x,y,'--og');
%The plot command again returns a handle, j

%TO DO 
%Change the color and width of the line using the 'set' command
set(j,'LineWidth',2)
set(j,'Color',[1 0 0])

%The objects in Matlab form a hierarchy.  The root object has one or more children which
%are figures. These figures have children which are sets of axes.  The axes have children
%which are lines.  You can retrieve the handle of the parent or child by extracting
%the 'parent' or 'children' field of the object using the 'get' command.

%TO DO 
%Retrieve the handle of the 
%axes that the line is plotted in by finding the parent of the line and
%storing it in a variable m.  You can see the properties by doing 'get(m)'
%Change the 'Box' property of the axis to 'off' (retain the quotes here)
m = get(j,'Parent')
set(m,'Box','off')
%TO DO 
%plot a new line of x vs z
%you should find that it completely replaced the previous curve
%to draw two curves on the same axis, you use the command 'hold on' before
%plotting the second curve.
hold on
j=plot(x,z)

%NOTE THAT you can get the handle to the current axis more conveniently,
%by using the command 'gca' - get current axis.  The command 'gcf' returns
%a handle to the current figure in the same way.

%=====================================================================

%Now let's investigate displaying images.

%TO DO 
%Create a 150 by 100 array of random numbers
%Display the array using the the command 'imagesc'
A = randn(150, 100)
imagesc(A)
%Use the command 'colorbar' to add a colorbar indicating what the colours
%mean
colorbar
%Use the command 'colormap(hot)' to change the colour map
colormap(hot)
%Use the command 'axis off' to remove the axis
axis off
%Use the command 'axis square' to make the pixels square
axis square
%TO DO 
%load in the image 'ucl.jpg' using the command 'imread'
%use the command 'whos' to note its dimensions and class
%display the image using the command 'whos'
A = imread('ucl.jpg')
whos
%display the image using the command 'image'
image

%Finally look up the following commands in the help
%'bar','surf','plot3','semilogx'







