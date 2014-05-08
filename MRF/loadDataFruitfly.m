
doEmbryo = 0;
doGut = 0;
doTouching = 1;

if doEmbryo
  X = imread('../data/embryo.png');
  X = imresize(X,0.25);
elseif doGut
  X = imread('../data/gut.png');
  X = imresize(X,0.5);
  X = histeq(X);
elseif doTouching
  X = imread('../data/touchingEmb.jpg');
  X = imresize(X,0.07);
end

X = mean(X,3);
X = (255 - X)/255;

figure;imagesc(X);
nNodes = length(X(:));
nStates = 2;
[nRows,nCols] = size(X);
