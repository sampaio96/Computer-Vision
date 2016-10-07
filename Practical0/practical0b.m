function r=practical0b

%Practical 0b. Math & debugging

%The goal of this practical is to investigate the math functions of
%matlab and to introduce the debugger

%Close all previous figures  
close all;
 

%Math in Matlab is basically set up to work with matrices. All common
%mathematical operations will work with arrays
A  = [1 2 ; 3 4]
B =  [2 0 ; 0 1]

%For example
C = A+B
D = sin(A)
E = exp(A)

%Note the difference between these commands.  What has happened in each
%case?
A
B
F = B*A
G = B.*A

%========================================================================

%The other important thing to note is that matlab math operations are MUCH 
%faster if they are done in parallel using array based commands than if we
%operate a for loop.

%You can time a set of commands by using the commands tic and toc
tic 
fprintf('Printing something to the screen\n');
toc

%TO DO 
%Create two arrays A and B of random numbers, each of size 400x400.  
%Time the following operations
%1. Creating array C by pointwise array division of A by B using the ./ operator
%2. Creating new array C by using two for loops to run through each location 
%   and dividing each element of A by each element of B in turn
tic
A = rand(400,400);
toc
tic
B = rand(400,400);
toc
tic
C = A./B;
toc
tic
    for a = 1:400
        for b = 1:400
            %C(a,b)=A(a,b)/B(a,b);
        end
    end
toc

%========================================================================

%You can stop a program anywhere by using the command 'keyboard' or by
%clicking next to the line number at the left.  Consider the code

A = randn(3,6);
A = A*6;
A = exp(-A);
A = A+2;

%TO DO 
%Stop the program between the first and second of these four lines and
%examining the contents of the matrix A

%TO DO 
%Step through the remaining lines by using the command "dbstep"

%Note - to continue the program running type dbcont
%     - to exit the debugger type dbquit
% Also see the keyboard shortcuts: F5, F10, F11.

% Note that you can select text in a .m file, and press F9 to execute just
% the hilighted command.

%Now move to practical0c






