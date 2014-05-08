
generateData;

%% Run MRF
% probA = ones(nNodes,nStates); %ones(1,nNodes)*(1/nStates);
% edgeStruct = createImgStructure(X,nStates);
% knownMuSigma = 0;
% if knownMuSigma  
%   [nodePot,edgePot] = producePotentials(X,edgeStruct,Mu,Sigma,probA);
%   edgeStruct.maxIter = max_iter; % Can set your own iterationmax
%    solveMRF;
% else
%   MuEst = [0.1,0.5];
%   SigmaEst = [0.1,0.1];
%   for i = 1:10
%     [nodePot,edgePot] = producePotentials(X,edgeStruct,MuEst,SigmaEst,probA);
%     [nodeBelMF,edgeBelMF,logZMF] = UGM_Infer_MeanField(nodePot,edgePot,edgeStruct);
%     [MuEst,SigmaEst] = em_max(nodeBelMF,X(:));
%     MuEst
%     SigmaEst
%   end
%   solveMRF
% end


%% Run Chan Vese
I = X;
% Customerized Mask
m = zeros(size(I,1),size(I,2));
m(10:size(I,1)-10,10:size(I,2)-10) = 1;
seg = chenvese(I,m,600,0.1); 

% Built-in Mask
seg = chenvese(I,'large',600,0.1); 
