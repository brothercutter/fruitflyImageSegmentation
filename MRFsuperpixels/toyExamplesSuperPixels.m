addpath(genpath('../UGM/'))
addpath(genpath('../MRF'))

%I = imread('../data/embryo.png');
%I = im2double(imresize(I,1));
%[nRows,nCols,dump] = size(I);
superpixel = 1; loadDataFruitfly; I = im2double(X);

doSuperPixels = false;
if doSuperPixels
  rng(215)
  createSuperPixels;
  numNbThreshold = 10;
  nStates = 2;
  createNeighbor;
  [edgePot,edgePotVec] = produceEdgePotentials(par.pb_emag,edgeStruct,idxArray,nStates);
  plotSuperPixelNeighbors;
  title('before merging small super pixels');
  mergeSmallSuperPixels;
  numNbThreshold = 10;
  createNeighbor;    
  nNodes = nSp;
  save([imgName,'sp.mat'],'Sp','nSp','nNodes','edgeStruct','idxArray','nStates','nRows','nCols','par','I_sp','numNbThreshold');
end
load([imgName,'sp.mat']);
[edgePot,edgePotVec] = produceEdgePotentials(par.pb_emag,edgeStruct,idxArray,nStates);
plotSuperPixelNeighbors;
title('after merging small super pixels');

% transform the pixel intensity to superpixel intensity
X = zeros(nSp,2);
for i = 1:nSp
  X(i,1) = 1-mean(I(find(Sp == i)));
  X(i,2) = std(I(find(Sp == i)));
end
X(:,1) = X(:,1)/max(X(:,1));

%{
%in_idx = [13,14,43,42,39,9,8,37,36,38,33,34,35,1];
in_idx = [17,21,16,22,20,10,33,4,31,8,34,9,7,32,27,6,27,1,26,45];
out_idx = setdiff(1:nSp,in_idx);

rng(210)
X(in_idx,1) = .5 + 0.1*randn(length(in_idx),1);
X(out_idx,1) = 0 + 0.1*randn(length(out_idx),1);;
%}

figure;
imagesc(labelSuperpixelsOnImage(Sp,X(:,1))); colorbar;

figure;
imagesc(labelSuperpixelsOnImage(Sp,X(:,2))); colorbar;

pause

if size(X,2) ==2
  MuEst = [0.5,0.8; 0.01,0.05];
  SigmaEst = [0.1,0.1; 0.01,0.01];
else 
  MuEst = [0,.3];
  SigmaEst = [.1,.1];
end

probA = ones(nNodes,nStates)*(1/nStates);
i = 0;
diff = 1;
iter_max = 20;


 while i<iter_max && diff>1e-5
     MuEst_o = MuEst;
     SigmaEst_o = SigmaEst;

     %[nodePot,edgePot] = producePotentials(X,edgeStruct,MuEst,SigmaEst,probA);
     nodePot = producePotentials(X,edgeStruct,MuEst,SigmaEst,probA);
     [nodeBelMF,edgeBelMF,logZMF] = UGM_Infer_MeanField(nodePot,edgePot,edgeStruct);
     %[nodeBelMF,edgeBelMF,logZMF] = UGM_Infer_LBP(nodePot,edgePot,edgeStruct);
     [MuEst,SigmaEst] = em_max(nodeBelMF,X);
     figure;
     imagesc(labelSuperpixelsOnImage(Sp,nodeBelMF(:,1)));
     title('Mean Field Estimates of Marginals');
     fprintf('(paused)\n');

     MuEst
     SigmaEst
     diff = max(norm(MuEst_o-MuEst),norm(SigmaEst_o - SigmaEst))
     pause
   end
   