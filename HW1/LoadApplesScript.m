% LoadApplesScript.m
% This optional script may help you get started with loading of photos and masks.
%
% Note that there are more elegant ways to organize your (photo, mask)
% pairs, and for sorting through many files in a directory. We don't use
% them here because we only have a small number of files, but consider
% useful functions like fileparts(). For simplicity, this example code just
% makes a few cell-arrays to hold the hard-coded filenames.

if( ~exist('apples', 'dir') || ~exist('testApples', 'dir') )
    display('Please change current directory to the parent folder of both apples/ and testApples/');
end

% Note that cells are accessed using curly-brackets {} instead of parentheses ().
n = 3; % Training+validation set size
Iapples = cell(n,1);
Iapples{1} = 'apples/Apples_by_kightp_Pat_Knight_flickr.jpg';
Iapples{2} = 'apples/ApplesAndPears_by_srqpix_ClydeRobinson.jpg';
Iapples{3} = 'apples/bobbing-for-apples.jpg';

IapplesMasks = cell(n,1);
IapplesMasks{1} = 'apples/Apples_by_kightp_Pat_Knight_flickr.png';
IapplesMasks{2} = 'apples/ApplesAndPears_by_srqpix_ClydeRobinson.png';
IapplesMasks{3} = 'apples/bobbing-for-apples.png';

Apples = cell(n-1,1);
AppleMasks = cell(n-1,1);
resizefactor = 0.2;
for iImage = 1:n-1
    tempApples = double(imresize(imread(Iapples{iImage}),resizefactor)) / 255;
    % Apples{iImage} is now a double-precision 3D matrix of size (width x height x 3). 
    % Each of the 3 color channels is now in the range [0.0, 1.0].
    
    Apples{iImage} = (reshape(tempApples,size(tempApples,1)*size(tempApples,2),3)).';
    % Apples{iImage} should be a matrix of size 3*(w*h)
    % size(Apples{iImage});
    
    tempAppleMasks = imresize(imread(IapplesMasks{iImage}),resizefactor);
    % These mask-images are often 3-channel, and contain grayscale values. We
    % would prefer 1-channel and just binary:
    tempAppleMasks = tempAppleMasks(:,:,2) > 128;  % Picked green-channel arbitrarily.
    AppleMasks{iImage} = (reshape(tempAppleMasks,size(tempAppleMasks,1)*size(tempAppleMasks,2),1)).';
    
    % Apples{iImage} should be a matrix of size 1*(w*h)    
    % size(AppleMasks{iImage});
end

Figure = horzcat(Apples{:}); % 3*I
Mask = horzcat(AppleMasks{:}); % 1*I

RGBApple = Figure( :, Mask );
RGBNonApple = Figure( :, ~Mask );

ValidationFig = imresize(imread(Iapples{n}),resizefactor);
ValidationMask = imresize(imread(IapplesMasks{n}),resizefactor);
tempVal = double(imresize(imread(Iapples{n}),resizefactor)) / 255;
imY = size(tempVal,1);
imX = size(tempVal,2);
FigureVal = (reshape(tempVal,imY*imX,3)).';
tempMaskVal = imresize(imread(IapplesMasks{n}),resizefactor);
tempMaskVal = tempMaskVal(:,:,2) > 128;
MaskVal = (reshape(tempMaskVal,size(tempMaskVal,1)*size(tempMaskVal,2),1)).';
RGBAppleVal = FigureVal( :, MaskVal );
RGBNonAppleVal = FigureVal( :, ~MaskVal );

% assert(isequal(length(Figure),length(RGBApple)+length(RGBNonApple)));
% assert(eq(size(RGBApple,2),sum(Mask)))
% assert(eq(size(RGBNonApple,2),size(Mask,2)-sum(Mask)))