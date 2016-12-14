function score = MeasurePatchSimilarityHere( Im2, pixelsTemplate, minY, minX )

SearchPadding = [0 0 0 0];
[cropCoordsFullI availPadW availPadN] = ...
    IndexesToSearchInFullImg(size(Im2), size(pixelsTemplate), [minY minX], SearchPadding);
subI = Im2(cropCoordsFullI(1):cropCoordsFullI(2), cropCoordsFullI(3):cropCoordsFullI(4), :);

% imgMaskCel = GetCelMaskPixels( maskSearchWhat, 1, minXs, maxXs, minYs, maxYs );

% Just use green channel for now (ie :,:,2).
[ncc movedSE scoreMaxs] = ...
    PatchDistPicker( pixelsTemplate(:,:,2), subI(:,:,2), ...
                     [availPadW availPadN], ...
                     'plain_nssd', ...   % builtinNCC, builtinNCCsearchInside, plain_ncc, maskedNCC, plain_nssd, masked_nssd
                     0 );
                 
if( numel(scoreMaxs) )
    score = max( 0, scoreMaxs(1) ); % Calling functions don't like negative scores.
else
    score = 0;
end;