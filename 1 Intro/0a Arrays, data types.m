function r=practical0a
%The goal of this practical is to introduce matlab programming
%It is divided into 3 sections

%Practical 0a. Arrays, data types

%Practical 0b. Math + the debugger

%Practical 0c. Data visualization

%You should work your way through this file, following the instructions
%marked "TO DO"

%===============================================================

%1. ARRAYS AND DATA TYPES

%close all previous figures
close all;

%clear the memory and all breakpoints
clear all;


%The most common data type in Matlab is an array of real numbers.
%To create an empty array, you can use the command zeros.  For example

A = zeros(3,4)

%creates an array full of zeros that has dimensions of 3 x 4.  Note that
%the first index refers to the vertical direction, not the horizontal
%direction.

%There are many other ways to initialize arrays that are full of ones,
%or random numbers

%TO DO  (YOU TYPE THE COMMANDS INTO THIS FILE)
%Initialize arrays B, C and D, of sizes (3x3), (4x6) and (2x3x2) using the
%commands 'zeros', 'rand' and 'randn' respectively.  To find out more about
%these commands type 'help zeros' etc. at the prompt

B = zeros(3,3)
C = rand(4,6)
D = rand(2,3,2)

%when you reach here, remove the command 'return' below 

%===============================================================

%We may not want the arrays to display to the command line every time:  if
%the array is very large, it will fill the screen up with numbers.  To
%suppress the screen output, append a semicolon to the end of the command.

%TO DO 
%create an array E of size (100x200) full of random numbers.  Suppress the
%output by appending a semicolon to the command that you used.

E = rand(100,200);

%===============================================================

%TO DO 
%Use the command 'whos' to investigate which commands are in memory

whos

%Notice that this command tells us about the "class" of the variable.  All
%of the arrays we have produced so far are of type "double array".
% Look at the workspace sub-window in the Matlab interface, and right-click
% there to see more about each variable.

%TO DO 
%Use the command 'size' to retreive the size of array A

size(A)
%===============================================================

%Let's create some other types of variable.

%String
F = 'Hello world'

%Unsigned integers (images will often come in this form)
G = zeros(10,10,'uint8')

%Cell arrays - a cell array is like an array of boxes.  Each box can hold 
%any other type of object - double array, string, or even another cell
%array. %To create a cell array:
H = cell(2,1);
%To index cell arrays (note curly-braces "{" and "}"  )
H{1,1} = 'Hi';
H{2,1} = zeros(4,2);

%TO DO 
%Use the command 'whos' to investigate which commands are in memory

whos
%===============================================================

%Now we'll investigate addressing arrays.  Lets clear all of the data
clear all;
%and create a single array
A = randn(5,7)
%There are 35 elements in this array.  We can index them in two different
%ways. 

%1. We can extract elements using full indexing using the format
B = A(2,3)

%TO DO 
%extract the bottom right element of the matrix A

C = A(5,7)

%We can set elements using a similar notation
A(1,5) = 10

%2. Alternatively, we can use linear indexing.  Here we use a single index
%to extract elements.  The elements are counted down the columns first, and
%then across the rows

B = A(4)
C = A(15)

%TO DO 
%Set the top right element of the matrix A to be equal to zero using single
%indexing.

A(31)=0

A

%===============================================================

%Now let's investigate some more complex ways to index arrays
%We can use the colon operator ':' to extract all of the elements in a
%given index position

%TO DO
%Investigate the result of the following commands (uncomment them in turn)
B = A(:,3)
C = A(2,:)
D = A(:)

%We can also extract just a subset of the array

%TO DO 
%Investigate the result of the following commands

B = A(2:4,:)
C = A(2:end,:)
D = A(2:3,4:6)


%We can also use the double colon notation to skip elements

%TO DO
%Investigate the result of the following commands

A

%B = A(1:3:end,:)
C = A(2,2:2:6)

%Finally, we can index one array with another
B = zeros(2,1);
B(1) = 4;
B(2) = 8

%TO DO
%Investigate the result of the following commands
C = A(B)


%===============================================================

%We can also reshape arrays to different sizes that have the same
%number of elements

A = randn(3,4);

%TO DO Use the command 'reshape' to change the shape of this array to be
%(2x6).  Note the order that the elements wind up in

reshape(A,2,6)


%===============================================================

%Finally, let's investigate concatenating arrays
A = randn(2,2)
B = randn(2,2)

%The square brackets notation allows us to concatenate arrays
%A space denotes horizontal concatenation, a semicolon denotes vertical 
%concatentation. 

%TO DO 
%Investigate the result of the following commands

C = [A B]
D = [A;B]
E = [1 2; 3 4]
F = [1 2:2:10; zeros(2,2) randn(2,4)]

%We can also concatenate arranys in arbitrary dimensions, using the command
%cat

%TO DO 
%Investigate the result of the following commands
C
G = cat(2,A,B)
H = cat(3,A,B)

%Finally, we can replicate arrays using the command repmat

%TO DO
%Investigate the result of the following commands
A
I = repmat(A,3,1)
B
J = repmat(B,2,3)


%Now move onto Practical 0b.
















