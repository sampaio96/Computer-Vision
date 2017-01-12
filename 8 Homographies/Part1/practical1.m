function practical1

    %The aim of practical 1 is to calculate the homography that best maps two
    %sets of points to one another.  We will (eventually) use this for creating
    %panoramas, and for calculating the 3d pose of planes.  You should use this
    %template for your code and fill in the missing sections marked "TO DO"

    %close all open figures
    close all;

    %set of two dimensional Cartesian points
    pts1Cart = [  240.5000   16.8351   33.5890  164.2696  149.1911;...
                  248.8770  193.5890   251.3901 168.4581  228.7723];

    %turn points to homogeneous representation
    pts1Hom = [pts1Cart; ones(1,size(pts1Cart,2))] % 3xI

    %define a homography
    H = [0.6 0.7 -100; 1.0 0.6 50; 0.001 0.002 1.0]; % 3x3
    %apply homography to points
    pts2Hom = H*pts1Hom

    %convert back to Cartesian
    pts2Cart = pts2Hom(1:2,:)./repmat(pts2Hom(3,:),2,1)

    %add a small amount of noise
    noiseLevel = 4.0;
    pts2Cart = pts2Cart+noiseLevel*randn(size(pts2Cart))

    %draw two set of two dimensional points
    %opens figure
    figure; set(gcf,'Color',[1 1 1]);
    %draw lines between each pair of points
    nPoint = size(pts1Cart,2);
    for (cPoint = 1:nPoint)
        %plot a green line between each pair of points
        plot([pts1Cart(1,cPoint) pts2Cart(1,cPoint)],[pts1Cart(2,cPoint) pts2Cart(2,cPoint)],'g-');
        %make sure we don't replace with next point
        hold on;
    end;

    %draws first set of points
    plot(pts1Cart(1,:),pts1Cart(2,:),'b.','MarkerSize',20);
    %remove axis
    set(gca,'Box','Off');

    %draws second set of points
    plot(pts2Cart(1,:),pts2Cart(2,:),'r.','MarkerSize',20);

    %now our goal is to transform the first points so that they map to the
    %second set of points
    
    HEst = calcBestHomography(pts1Cart, pts2Cart);

    %now we will see how well the routine works by applying the mapping and
    %measuring the square  distance between the desired and actual positions

    %apply homography to points
    pts2EstHom = HEst*pts1Hom

    %convert back to Cartesian
    pts2EstCart = pts2EstHom(1:2,:)./repmat(pts2EstHom(3,:),2,1)

    %calculate mean squared distance from actual points
    sqDiff = mean(sum((pts2Cart-pts2EstCart).^2))

    %draw figure with points before and after
    %draw two set of two dimensional points
    %opens figure
    figure; set(gcf,'Color',[1 1 1]);
    %draw lines between each pair of points
    nPoint = size(pts1Cart,2)
    for (cPoint = 1:nPoint)
        %plot a green line pairs of actual and estimated points
        plot([pts2Cart(1,cPoint) pts2EstCart(1,cPoint)],[pts2Cart(2,cPoint) pts2EstCart(2,cPoint)],'g-');
        %make sure we don't replace with next point
        hold on;
    end;

    %draws second set of points
    plot(pts2Cart(1,:),pts2Cart(2,:),'r.','MarkerSize',20);
    %remove axis
    set(gca,'Box','Off');

    %draws estimated positions of second set of points
    plot(pts2EstCart(1,:),pts2EstCart(2,:),'m.','MarkerSize',20);

    %other things **** TO DO ****
    %1. Convince yourself that the homography is ambiguous up to scale (by
    %multiplying it by a constant factor and showing it does the same thing).
    %Can you see why this is mathematically the case?
    % Yes, because x is obtained from u and v by multiplying and dividing by
    % the same number of parameters in hom matrix.
    %2. Show empirically that your homography routine can EXACTLY map any four points to any
    %other four points
    %3. Now move to practical 1b.
    
end

%==========================================================================
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
    h = solveAXEqualsZero(A);

    %reshape h into the matrix H
    H = reshape(h,3,3).';
end

%==========================================================================
function x = solveAXEqualsZero(A)
    [U,L,V] = svd(A);
    x = V(:, end);
end