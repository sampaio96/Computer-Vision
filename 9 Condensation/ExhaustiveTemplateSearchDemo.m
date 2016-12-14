Im1 = double( imread( 'HillsRdSkipFrames_0000871.png' ));
Im2 = double( imread( 'HillsRdSkipFrames_0000878.png' ));
%imagesc(Im1/255);

% User will need to click-drag the corners of a rectangle on the current
% figure. Hit <enter> when done, so that coordinates are saved in pos:
% [upper left's x, y, width, height]
%H = imrect;
%pos = wait(H);
% OR: 
pos = [253.81       256.07       17.419       25.714 ];
pos = round(pos);


% These steps are QUITE round-about for isolating the "template" pixels we
% wish to track, but makes it easy to see how one could select a more
% freehand shape than a rectangle, like an impoly, imellipse, or
% imfreehand.
maskSearchWhat = uint8(zeros(size(Im1, 1), size(Im1, 2), 1));
maskSearchWhat( pos(2):pos(2)+pos(4), pos(1):pos(1)+pos(3) ) = 1;
clear H pos;
[Ys Xs] = find(maskSearchWhat == 1);
% Now we have all the coord's of this blob.

% Two equivalent ways to grab the sub-image:
minXs = min(Xs); 
maxXs = max(Xs);    
minYs = min(Ys); 
maxYs = max(Ys);
pixelsTemplate = Im1(minYs:maxYs,   minXs:maxXs, :);  % Need colRange, rowRange
% OR
% upper-left corner, width, height:
% rect_Cel = [min(Xs)     min(Ys)     max(Xs)-min(Xs)    max(Ys)-min(Ys)];
% pixelsTemplate = imcrop( I, rect_Cel );  % NOTE: imcrop uses (x,y) NOT (col, row)
% imagesc(pixelsTemplate/255)

contour = bwtraceboundary(logical(maskSearchWhat), [Ys(1), Xs(1)], 'NE', 8, inf, 'counterclockwise');
% figure; imagesc(Im1/255);
% hold on;
% plot(contour(:,2),contour(:,1),'g','LineWidth',2);

% Saving these off for the Lab practical
% minY = minYs; minX = minXs;
% save('Template', 'pixelsTemplate', 'minY', 'minX')

SearchPadding = [50 50 50 50];
[cropCoordsFullI availPadW availPadN] = ...
    IndexesToSearchInFullImg(size(Im2), size(pixelsTemplate), [minYs minXs], SearchPadding);
subI = Im2(cropCoordsFullI(1):cropCoordsFullI(2), cropCoordsFullI(3):cropCoordsFullI(4), :);

imgMaskCel = GetCelMaskPixels( maskSearchWhat, 1, minXs, maxXs, minYs, maxYs );

% Just use green channel for now (ie :,:,2).
[ncc movedSE scoreMaxs] = ...
    PatchDistPicker( pixelsTemplate(:,:,2), subI(:,:,2), ...
                     [availPadW availPadN], ...
                     'builtinNCCsearchInside', ...   % builtinNCC, builtinNCCsearchInside, plain_ncc, maskedNCC, plain_nssd, masked_nssd
                     imgMaskCel );

                 
% Visualizations
% ================================

%     imagesc(ncc)
%     colormap gray
%     pause(2)
%     surf(ncc), shading flat

hImg = figure;
set(0,'CurrentFigure',hImg)
set(gcf,'Position',[20 300 481 361]);
set(gcf,'Color',[1 1 1]);
set(gca,'Box','Off');
title('Locations');

% Put down the 2nd (i.e. current frame)
imagesc(Im2/255);
hold on;

% Plot a green box showing where the highest-likelihood match was found
bestSouth = movedSE(1,1);
bestEast = movedSE(1,2);
plot(contour(:,2) + bestEast, contour(:,1) + bestSouth,'g','LineWidth',2);

% Plot in yellow: where the search was centered
s = size(pixelsTemplate);
ul = [minYs, minXs];
ll = [minYs + s(1), minXs];
lr = [minYs + s(1), minXs + s(2)];
ur = [minYs, minXs + s(2)];
boundingBox = [ul; ll; lr; ur; ul];
plot( boundingBox(:,2), boundingBox(:,1), 'y:', 'LineWidth', 1 );

% Plot in blue: the outer limits of where the likelihood was computed
ul = [cropCoordsFullI(1), cropCoordsFullI(3)];
ll = [cropCoordsFullI(2), cropCoordsFullI(3)];
lr = [cropCoordsFullI(2), cropCoordsFullI(4)];
ur = [cropCoordsFullI(1), cropCoordsFullI(4)];
boundingBox = [ul; ll; lr; ur; ul];
plot( boundingBox(:,2), boundingBox(:,1), 'b:', 'LineWidth', 3 );

hold off

