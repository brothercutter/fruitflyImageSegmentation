clear all
close all
addpath(genpath('../'));
%generateData;
loadDataFruitfly;

loadMRI = 0;
if loadMRI
  X = imread('../Chan-Vese/brain.jpg');
  X = mean(X,3);
  X = X/255 + 0.01*randn(size(X,1),size(X,2));
  figure;imagesc(X);colorbar 
  nNodes = length(X(:));
  nStates = 2;
  [nRows,nCols] = size(X);
end

 MuEst = [0.5,1];
 SigmaEst = [0.1,0.3];
 i = 0;
 diff = 1;
 iter_max = 20;

%% Run MRF

%  probA = ones(nNodes,nStates)*(1/nStates);
%  edgeStruct = createImgStructure(X,nStates);
% 
% 
% 
%  while i<iter_max && diff>1e-3
%      MuEst_o = MuEst;
%      SigmaEst_o = SigmaEst;
% 
%      [nodePot,edgePot] = producePotentials(X,edgeStruct,MuEst,SigmaEst,probA);
%      [nodeBelMF,edgeBelMF,logZMF] = UGM_Infer_MeanField(nodePot,edgePot,edgeStruct);
%      [MuEst,SigmaEst] = em_max(nodeBelMF,X(:));
%      figure;
%      imagesc(reshape(nodeBelMF(:,2),nRows,nCols));
%      colormap gray
%      title('Mean Field Estimates of Marginals');
%      fprintf('(paused)\n');
%      
%      MuEst
%      SigmaEst
%      diff = max(norm(MuEst_o-MuEst),norm(SigmaEst_o - SigmaEst))
%    end
%    solveMRF
% 

% stop

%% Run Chan Vese

I = X;
% Customerized Mask
m = zeros(size(I,1),size(I,2));
m(10:size(I,1)-10,10:size(I,2)-10) = 1;
seg = chanvese(I,m,600,0.1,[],0,[],0,1e-5); 

Y0 = X(find(seg == 0));
Y1 = X(find(seg == 1));

MuEst(1) = mean(Y0);
MuEst(2) = mean(Y1);

SigmaEst(1) = std(Y0);
SigmaEst(2) = std(Y1);


% Built-in Mask
% seg = chanvese(I,'large',1000,0.1,[],0); 
