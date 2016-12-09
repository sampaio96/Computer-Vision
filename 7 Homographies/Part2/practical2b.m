function practical2b

%The goal of this part of the practical is to take a real image containing
%a planar black square and figure out the transformation between the square
%and the camera.  We will then draw a wire-frame cube with it's base
%corners at the corner of the square.  You should use this
%template for your code and fill in the missing sections marked "TO DO"

%load in image 
im = imread('test104.jpg');

%define points on image
xImCart = [  140.3464  212.1129  346.3065  298.1344   247.9962;...
             308.9825  236.7646  255.4416  340.7335   281.5895];
         
%define 3D points of plane
XCart = [-50 -50  50  50 0 ;...
          50 -50 -50  50 0;...
           0   0   0   0 0];

%We assume that the intrinsic camera matrix K is known and has values
K = [640  0    320;...
     0    640  240;
     0    0    1];

%draw image and 2d points
figure; set(gcf,'Color',[1 1 1]);
imshow(im); axis off; axis image; hold on;
plot(xImCart(1,:),xImCart(2,:),'r.','MarkerSize',10);
       
%TO DO Use your routine to calculate TEst, the extrinsic matrix relating the
%plane position to the camera position.

TEst = estimatePlanePose(xImCart,XCart,K);

%define 3D points of plane
XWireFrameCart = [-50 -50  50  50 -50 -50  50  50;...
                   50 -50 -50  50  50 -50 -50  50;...
                    0   0   0   0 -100 -100 -100 -100];

%TO DO Draw a wire frame cube, by projecting the vertices of a 3D cube
%through the projective camera and drawing lines betweeen the resulting 2d image
%points

CV = projectiveCamera(K,TEst,XWireFrameCart);
% A = CubeVertices(1, 2; 1, 4; 1, 5; 2, 3; 2, 6; 3, 4; 3, 7; 5, 6; 5, 8; 6, 7; 7, 8)
x = [CV(1,1), CV(1,2); CV(1,1), CV(1,4); CV(1,1), CV(1,5); CV(1,2), CV(1,3); CV(1,2), CV(1,6); CV(1,3), CV(1,4); CV(1,3), CV(1,7); CV(1,5), CV(1,6); CV(1,5), CV(1,8); CV(1,6), CV(1,7); CV(1,7), CV(1,8)].';
y = [CV(2,1), CV(2,2); CV(2,1), CV(2,4); CV(2,1), CV(2,5); CV(2,2), CV(2,3); CV(2,2), CV(2,6); CV(2,3), CV(2,4); CV(2,3), CV(2,7); CV(2,5), CV(2,6); CV(2,5), CV(2,8); CV(2,6), CV(2,7); CV(2,7), CV(2,8)].';
plot(x,y,'-or');

%QUESTIONS TO THINK ABOUT...

%Do the results look realistic?
%If not, then what factors do you think might be causing this?
% No shadows/filling, and picture is curved (panoramic), while the cube
% isn't. The top plane is also colored black.

function xImCart = projectiveCamera(K,T,XCart);

%replace this
%xImCart = [];

%TO DO convert Cartesian 3d points XCart to homogeneous coordinates XHom
XHom=[XCart;ones(1,size(XCart,2))];
%TO DO apply extrinsic matrix to XHom to move to frame of reference of
%camera
TX=T*XHom;
%TO DO project points into normalized camera coordinates xCamHom by (achieved by
%removing fourth row)
xCamHom=TX(1:end-1,:);
%TO DO move points to image coordinates xImHom by applying intrinsic matrix
xImHom=K*xCamHom;
%TO DO convert points back to Cartesian coordinates xImCart
xImCart = xImHom(1:2,:)./repmat(xImHom(3,:),2,1);

function T = estimatePlanePose(xImCart,XCart,K)

%replace this
%T = [];

%TO DO Convert Cartesian image points xImCart to homogeneous representation
xImHom=[xImCart;ones(1,size(xImCart,2))];

%TO DO Convert image co-ordinates xImHom to normalized camera coordinates
xCamHom=K\xImHom;

%TO DO Estimate homography H mapping homogeneous (x,y)
%coordinates of positions in real world to xCamHom.  Use the routine you wrote for
%Practical 1B.
H = calcBestHomography(XCart, xCamHom)
%TO DO Estimate first two columns of rotation matrix R from the first two
%columns of H using the SVD
[U L V]=svd(H(:,1:2));
O=U*[1 0; 0 1; 0 0]*V';
R=O(:,1:2);
%TO DO Estimate the third column of the rotation matrix by taking the cross
%product of the first two columns
R=[R cross(R(:,1),R(:,2))];
%TO DO Check that the determinant of the rotation matrix is positive - if
%not then multiply last column by -1.
if det(R) < 0
    R = R*[1 0 0; 0 1 0; 0 0 -1];
end
%TO DO Estimate the translation t by finding the appropriate scaling factor k
%and applying it to the third colulmn of H
k = sum(sum(R(:,1:2)./H(:,1:2)))/6;
t = k*H(:,3);

%TO DO Check whether t_z is negative - if it is then multiply t by -1 and
%the first two columns of R by -1.
if t(3)<0
    t = -t;
    R = R*[-1 0 0; 0 -1 0; 0 0 -1];
end
%assemble transformation into matrix form
T  = [R t;0 0 0 1];

function H = calcBestHomography(pts1Cart, pts2Cart)

%should apply direct linear transform (DLT) algorithm to calculate best
%homography that maps the points in pts1Cart to their corresonding matchin in 
%pts2Cart

%****TO DO ****: replace this
%H = eye(3);

%**** TO DO ****;
%first turn points to homogeneous
% pts1Hom = [pts1Cart; ones(1,size(pts1Cart,2))];
% pts2Hom = [pts2Cart; ones(1,size(pts2Cart,2))];

%then construct A matrix which should be (10 x 9) in size
A=[];
for i=1:size(pts1Cart,2)
    u = pts1Cart(1,i);
    v = pts1Cart(2,i);
    x = pts2Cart(1,i);
    y = pts2Cart(2,i);
    A=[A;
        0 0 0 -u -v -1 y*u y*v y;
        u v 1 0 0 0 -x*u -x*v -x];
end

%solve Ah = 0 by calling
h = solveAXEqualsZero(A)
%(you have to write this routine too - see below)
%reshape h into the matrix H
H=reshape(h,3,3)';
%Beware - when you reshape the (9x1) vector x to the (3x3) shape of a homography, you must make
%sure that it is reshaped with the values going first into the rows.  This
%is not the way that the matlab command reshape works - it goes columns
%first.  In order to resolve this, you can reshape and then take the
%transpose


%==========================================================================
function x = solveAXEqualsZero(A);

%****TO DO **** Write this routine 
[U,S,V] = svd(A);
x=V(:,end);