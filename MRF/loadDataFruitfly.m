X = imread('../data/embryo.png');
X = imresize(X,0.25);
X = mean(X,3);
X = (255 - X)/255;
figure;imagesc(X);
nNodes = length(X(:));
nStates = 2;
[nRows,nCols] = size(X);
