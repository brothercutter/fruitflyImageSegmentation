addpath(genpath('../'));
%generateData;
%loadDataFruitfly;

loadMRI = 1;
if loadMRI
  X = imread('../Chan-Vese/brain.jpg');
  X = mean(X,3);
  X = X/255 + 0.01*randn(size(X,1),size(X,2));
  figure;imagesc(X);colorbar 
  nNodes = length(X(:));
  nStates = 2;
  [nRows,nCols] = size(X);
end


%% Run MRF

 probA = ones(nNodes,nStates)*(1/nStates);
 edgeStruct = createImgStructure(X,nStates);

 MuEst = [0,1];
 SigmaEst = [0.1,0.3];
   for i = 1:20
     [nodePot,edgePot] = producePotentials(X,edgeStruct,MuEst,SigmaEst,probA);
     [nodeBelMF,edgeBelMF,logZMF] = UGM_Infer_MeanField(nodePot,edgePot,edgeStruct);
     [MuEst,SigmaEst] = em_max(nodeBelMF,X(:));
     MuEst

     SigmaEst
     figure;
     imagesc(reshape(nodeBelMF(:,2),nRows,nCols));
     colormap gray
     title('Mean Field Estimates of Marginals');
     fprintf('(paused)\n');

   end
   solveMRF


% stop

%% Run Chan Vese

I = X;
I = histeq(X);
% Customerized Mask
m = zeros(size(I,1),size(I,2));
m(10:size(I,1)-10,10:size(I,2)-10) = 1;
seg = chanvese(I,m,600,0.1); 

% Built-in Mask
seg = chenvese(I,'large',600,0.1); 
