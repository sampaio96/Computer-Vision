function practical1B

%the aim of the second part of practical 1 is to use the homography routine
%that you established in the first part of the practical.  We are going to
%make a panorama of several images that are related by a homography.  I
%provide 3 images (one of which is has a large surrounding region) and a
%matching set of points between these images.

%close all open figures
close all;

%load in the required data
load('PracticalData','im1','im2','im3','pts1','pts2','pts3','pts1b');
%im1 is center image with grey background
%im2 is left image 
%pts1 and pts2 are matching points between image1 and image2
%im3 is right image
%pts1b and pts3 are matching points between image 1 and image 3
im1b = im1;
%show images and points
figure; set(gcf,'Color',[1 1 1]);image(uint8(im1));axis off;hold on;axis image;
plot(pts1(1,:),pts1(2,:),'r.'); 
plot(pts1b(1,:),pts1b(2,:),'m.');
figure; set(gcf,'Color',[1 1 1]);image(uint8(im2));axis off;hold on;axis image;
plot(pts2(1,:),pts2(2,:),'r.'); 
figure; set(gcf,'Color',[1 1 1]);image(uint8(im3));axis off;hold on;axis image;
plot(pts3(1,:),pts3(2,:),'m.'); 


H12 = calcBestHomography(pts1,pts2);
H13 = calcBestHomography(pts1b,pts3);

%==========================================================================

%for every pixel in image 2
for y = 1:size(im1,1)
    for x = 1:size(im1,2)
    
    %transform this pixel position to find where it is in the coordinates of image 1 
        A = (H12)*[x;y;1];
        pixel1 = round([A(1)/A(3);A(2)/A(3)]);
        if pixel1(1) <= size(im2,2) && pixel1(1) > 0 && pixel1(2) <= size(im2,1) && pixel1(2) > 0
            %copy pixel colour from image 2 pixel to current position in image 1 
            im1(y,x,:) = im2(pixel1(2),pixel1(1),:);
        end
        
        A = (H13)*[x;y;1];
        pixel1 = round([A(1)/A(3);A(2)/A(3)]);
        if pixel1(1) <= size(im2,2) && pixel1(1) > 0 && pixel1(2) <= size(im2,1) && pixel1(2) > 0
            %copy pixel colour from image 2 pixel to current position in image 1 
            im1(y,x,:) = im3(pixel1(2),pixel1(1),:);
        end
        
    end
end

figure; set(gcf,'Color',[1 1 1]);image(uint8(im1));axis off;hold on;axis image;

function H = calcBestHomography(pts1Cart, pts2Cart)

    %should apply direct linear transform (DLT) algorithm to calculate best
    %homography that maps the points in pts1Cart to their corresonding matching in 
    %pts2Cart
    H = eye(3);
    nPoint = size(pts1Cart,2);

    %first turn points to homogeneous
    pts1 = [pts1Cart; ones(1,nPoint)]; %3x5
    pts2 = [pts2Cart; ones(1,nPoint)]; %3x5
    
    %then construct A matrix which should be (10 x 9) in size
    A = zeros(nPoint*2,9);
    for i = 1:nPoint
        A(2*i-1:2*i,:) = [0, 0, 0, -pts1(:,i)', pts2(2,i)*pts1(:,i)'; pts1(:,i)', 0, 0, 0, -pts2(1,i)*pts1(:,i)'];
    end

    %solve Ah = 0 by calling
    [U,L,V] = svd(A);
    h = V(:, end);

    %reshape h into the matrix H
    H = reshape(h,3,3).';