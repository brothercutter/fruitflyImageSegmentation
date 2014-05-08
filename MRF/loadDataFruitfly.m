doEmbryo = 0;
doGut = 1;
if doEmbryo
  X = imread('../data/embryo.png');
  X = imresize(X,0.25);
elseif doGut
  X = imread('../data/gut.png');
  X = imresize(X,0.5);
end

X = mean(X,3);
X = (255 - X)/255;
X = histeq(X);
figure;imagesc(X);
nNodes = length(X(:));
nStates = 2;
[nRows,nCols] = size(X);
