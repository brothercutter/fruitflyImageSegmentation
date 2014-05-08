
%generateData;
loadDataFruitfly;


%% Run MRF

 probA = ones(nNodes,nStates)*(1/nStates);
 edgeStruct = createImgStructure(X,nStates);

 MuEst = [0.1,0.5];
 SigmaEst = [0.1,0.3];
   for i = 1:50
     [nodePot,edgePot] = producePotentials(X,edgeStruct,MuEst,SigmaEst,probA);
     [nodeBelMF,edgeBelMF,logZMF] = UGM_Infer_MeanField(nodePot,edgePot,edgeStruct);
     [MuEst,SigmaEst] = em_max(nodeBelMF,X(:));
     MuEst
     SigmaEst
 
   end
   solveMRF


% stop

%% Run Chan Vese
%I = X;
I = histeq(X);
% Customerized Mask
m = zeros(size(I,1),size(I,2));
m(10:size(I,1)-10,10:size(I,2)-10) = 1;
seg = chanvese(I,m,600,0.1); 

% Built-in Mask
seg = chenvese(I,'large',600,0.1); 
