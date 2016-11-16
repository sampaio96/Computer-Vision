function r=practicalMixGaussA

%load in test image and ground truth.
resizefactor = 0.5;
im = imresize(imread('apples/Apples_by_kightp_Pat_Knight_flickr.jpg'),resizefactor);
gt = imresize(imread('apples/Apples_by_kightp_Pat_Knight_flickr.png'),resizefactor);

% Load training data RGBApple, RGBNonApple, 3*I vectors
run('LoadApplesScript.m');

%display test image and ground truth;
close all;
figure; set(gcf,'Color',[1 1 1]);
subplot(1,3,1); imagesc(im); axis off; axis image;
subplot(1,3,2); imagesc(gt); colormap(gray); axis off; axis image;
drawnow;

%fit Gaussian model for apple data
mixGaussApple = fitMixGauss(RGBApple(:,randi([1 size(RGBApple,2)],1,10000)),10);

%fit Gaussian model for non-apple data
mixGaussNonApple = fitMixGauss(RGBNonApple(:,randi([1 size(RGBNonApple,2)],1,10000)),20);

[imY imX imZ] = size(im);
data = double(reshape(im,imY*imX,imZ).') / 255;

%determine posterior for each pixel being apple
[imY imX imZ] = size(im);

posteriorApple = zeros(imY,imX);
for (cY = 1:imY); 
    fprintf('Processing Row %d\n',cY);     
    for (cX = 1:imX);          
        %extract this pixel data
        thisPixelData = squeeze(double(im(cY,cX,:)));
        %calculate likelihood of this data given apple model
        likeApple = calcGaussianProb(thisPixelData,mixGaussApple);
        %calculate likelihood of this data given non-apple model
        likeNonApple = calcGaussianProb(thisPixelData,mixGaussNonApple);
        %calculate posterior probability
        num = likeApple*0.3;
        den = num + likeNonApple*0.7;
        posteriorApple(cY,cX) = num/den;
    end;
end;

%draw skin posterior
clims = [0, 1];
subplot(1,3,3); imagesc(posteriorApple, clims); colormap(gray); axis off; axis image;
% set(gca, 'clim', [0, 1]);

%==========================================================================
%==========================================================================

function bound = getBound(data,mixGaussEst,responsibilities)

%find total number of data items
nData = size(data,2);

%initialize bound
bound = 0;

%run through each data item
for(cData = 1:nData)
    %extract this data
    thisData = data(:,cData);    
    %extract this q(h)
    thisQ = responsibilities(:,cData);
    
    nDim = length(thisData);
    assert(eq(nDim,mixGaussEst.d));
    
    %TO DO - calculate contribution to bound of this datapoint
    %Replace this
    boundValue = 0;    
    for(cGauss = 1:mixGaussEst.k)
        boundValue = boundValue + thisQ(cGauss) * log ( (mixGaussEst.weight(cGauss) * 1/((2*pi)^(nDim/2)*(norm(mixGaussEst.cov(:,:,cGauss)))^(1/2)) * exp(-0.5*((thisData-mixGaussEst.mean(cGauss)).')*inv(mixGaussEst.cov(:,:,cGauss))*(thisData-mixGaussEst.mean(cGauss)))) / thisQ(cGauss));
    end
    
    %add to total log like
    bound = bound+boundValue;
end;

function logLike = getMixGaussLogLike(data,mixGaussEst)

[nDim nData] = size(data);

%initialize log likelihoods
logLike = 0;

xx = reshape(data,nDim,1,nData);
% assert(isequal(size(xx),[nDim 1 nData]));
xx = repmat(xx,1,mixGaussEst.k,1);
% assert(isequal(size(xx),[nDim mixGaussEst.k nData]));
xx = xx - mixGaussEst.mean;
% assert(isequal(size(xx),[nDim mixGaussEst.k nData]));
xx = mat2cell(xx, nDim, mixGaussEst.k, ones(1,nData));

ol = cellfun(@ (A) (arrayfun((@ (cGauss) (A(:,cGauss)' * mixGaussEst.cov(:,:,cGauss) * A(:,cGauss))), (1:mixGaussEst.k))),xx,'un',0);
% assert(isequal(size(ol{1}),[1 mixGaussEst.k]));
ol = cellfun(@ (A) (-0.5 * exp(A)), ol,'un',0);
% assert(isequal(size(ol{1}),[1 mixGaussEst.k]));
oj = arrayfun((@ (cGauss) (mixGaussEst.weight(1) / ((2*pi)^(nDim/2) * norm(mixGaussEst.cov(1))^(1/2)))), (1:mixGaussEst.k));
ol = cellfun(@ (A) (oj .* A), ol,'un',0);
% assert(isequal(size(ol{1}),[1 mixGaussEst.k]));
oi = cellfun(@ (A) (A(:,1) + A(:,2)), ol,'un',0);
% assert(isequal(size(oi{1}),[1 1]));
    
% like = sum(cell2mat(oi));
oi = cellfun(@ (A) (log(A)), oi,'un',0);
logLike = sum(cell2mat(oi)); %scalar

function like = calcGaussianProb(data, mixGauss)

D = length(data);
for cGauss = 1:mixGauss.k
    like = mixGauss.weight(cGauss) * (1/((2*pi)^(D/2)*(norm(mixGauss.cov(cGauss)))^(1/2)))*exp(-0.5*((data-mixGauss.mean(cGauss)).')*inv(mixGauss.cov(cGauss))*(data-mixGauss.mean(cGauss)));
end

% function [meanData covData] = fitGaussianModel(data)
% 
% [nDim nData] = size(data);
% 
% %TO DO (a): replace this
% meanData = mean(data.').'; % mean of each column gives out a 3D mean
% 
% xs = data - meanData;
% covData = xs * xs';
% covData = covData/(nData-1);

% assert(isequal(size(meanData),[3 1]));
% assert(isequal(size(covData),[3 3]));
% assert(isequal(covData,cov(data.'))); % is fine, but doesn't need to be
% computed each time

function mixGaussEst = fitMixGauss(data,k)

[nDim nData] = size(data);
mixGaussEst.d = nDim;
mixGaussEst.k = k;
mixGaussEst.weight = (1/k)*ones(1,k);
mixGaussEst.mean = 2*randn(nDim,k);
for (cGauss =1:k)
    mixGaussEst.cov(:,:,cGauss) = (0.5+1.5*rand(1))*eye(nDim,nDim);
end;

%calculate current likelihood and bound
logLike = getMixGaussLogLike(data,mixGaussEst);
fprintf('Log Likelihood Iter 0 : %4.3f\n',logLike);

nIter = 5; % updates weights, and means and cov too to increase likelihood
for (cIter = 1:nIter)

% Expectation step
xx = reshape(data,3,1,nData);
% assert(isequal(size(xx),[nDim 1 nData]));
xx = repmat(xx,1,k,1);
% assert(isequal(size(xx),[nDim k nData]));
xx = xx - mixGaussEst.mean;
% assert(isequal(size(mixGaussEst.mean),[nDim k]));  
% assert(isequal(size(xx),[nDim k nData]));
xx = mat2cell(xx, nDim, k, ones(1,nData));
% assert(isequal(size(xx{1}),[nDim k]));

ol = cellfun(@ (A) (arrayfun((@ (cGauss) (A(:,cGauss)' * mixGaussEst.cov(:,:,cGauss) * A(:,cGauss))), (1:mixGaussEst.k))),xx,'un',0);
% assert(isequal(size(ol{1}),[1 k]));
ol = cellfun(@ (A) (-0.5 * exp(A)), ol,'un',0);
% assert(isequal(size(ol{1}),[1 k]));
oj = arrayfun((@ (cGauss) (mixGaussEst.weight(1) / ((2*pi)^(nDim/2) * norm(mixGaussEst.cov(1))^(1/2)))), (1:mixGaussEst.k));
% assert(isequal(size(oj),[1 k]));
ol = cellfun(@ (A) (oj .* A), ol,'un',0);
% assert(isequal(size(ol{1}),[1 k]));

oi = reshape(cell2mat(ol),mixGaussEst.k,nData);
% assert(isequal(size(oi),[k nData]));
postHidden = oi ./ sum(oi);
% assert(isequal(size(postHidden),[k nData]))

% Maximization Step

%for each constituent Gaussian
sumr = sum(sum(postHidden(:,:)));
for (cGauss = 1:k) 
    sumIr = sum(postHidden(cGauss,:));

    mixGaussEst.weight(cGauss) = sumIr/sumr;

    ji = data(:,:) .* postHidden(cGauss,:);
    mixGaussEst.mean(:,cGauss) = sum(ji,2) ./ sumIr;

    jj = data - mixGaussEst.mean(:,cGauss);
    % assert(isequal(size(jj),[nDim 1]));
    jk = mat2cell(reshape(jj,mixGaussEst.d,1,nData), nDim, 1, ones(1,nData));
    jk = cellfun(@ (A) (A * A'), jk,'un',0);
    jk = reshape(cell2mat(jk),nDim,nDim,nData) .* reshape(postHidden(cGauss,:),1,1,nData);
    mixGaussEst.cov(:,:,cGauss) = sum(jk,3) ./ sumIr;
end;

%calculate the log likelihood
logLike = getMixGaussLogLike(data,mixGaussEst);
fprintf('Log Likelihood Iter %d : %4.3f\n',cIter,logLike);

% bound = getBound(data,mixGaussEst,postHidden);
% fprintf('Bound Iter %d : %4.3f\n',cIter,bound);
end;