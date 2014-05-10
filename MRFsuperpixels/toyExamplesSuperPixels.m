addpath(genpath('../UGM/'))
addpath(genpath('../MRF'))

I = imread('../data/embryo.png');
I = im2double(imresize(I,1));
[nRows,nCols] = size(I);

%createSuperPixels;
nSp = length(unique(Sp(:)));

% transform the pixel intensity to superpixel intensity
X = zeros(nSp,1);
for i = 1:nSp
  X(i) = mean(I(find(Sp == i)));
end

X = 1 - X;
figure;
imagesc(labelSuperpixelsOnImage(Sp,X));
nStates = 2;
nNodes = nSp;
createNeighbor;

MuEst = [0.1,.9];
SigmaEst = [0.2,0.1];
probA = ones(nNodes,nStates)*(1/nStates);
i = 0;
diff = 1;
iter_max = 20;


 while i<iter_max && diff>1e-5
     MuEst_o = MuEst;
     SigmaEst_o = SigmaEst;

     [nodePot,edgePot] = producePotentials(X(:),edgeStruct,MuEst,SigmaEst,probA);
     [nodeBelMF,edgeBelMF,logZMF] = UGM_Infer_MeanField(nodePot,edgePot,edgeStruct);
     [MuEst,SigmaEst] = em_max(nodeBelMF,X(:));
     figure;
     imagesc(labelSuperpixelsOnImage(Sp,nodeBelMF(:,1)));
     colormap gray
     title('Mean Field Estimates of Marginals');
     fprintf('(paused)\n');

     MuEst
     SigmaEst
     diff = max(norm(MuEst_o-MuEst),norm(SigmaEst_o - SigmaEst))
     pause
   end
   