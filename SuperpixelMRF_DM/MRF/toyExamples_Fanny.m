clear all

close all

%% Combined MRF and Level Set
%max_iter = 300;
%generateData;
superpixel = 0;
loadDataFruitfly;
Epsilon = 1;

iter_max = 50;
mf_iter_max = 500;
i = 0;
diff = 1;

% %% Run MRF
% probA = ones(nNodes,nStates); %ones(1,nNodes)*(1/nStates);
% edgeStruct = createImgStructure(X,nStates);
% knownMuSigma = 0;
% if knownMuSigma  
%    [nodePot,edgePot] = producePotentials(X,edgeStruct,Mu,Sigma,probA);
%    edgeStruct.maxIter = max_iter; % Can set your own iterationmax
%    solveMRF;
% else
%   MuEst = [0.5 0.99];  % For Flowers Embryo, 0.1, 0.5 worked ok. % brain: 0.7,0.9
%   SigmaEst = [0.001,0.001]; % brain: 0.3 0.01
%   while i<iter_max && diff>1e-3
%     MuEst_o = MuEst;
%     SigmaEst_o = SigmaEst;
%     [nodePot,edgePot] = producePotentials(X,edgeStruct,MuEst,SigmaEst,probA);
%     [nodeBelMF,edgeBelMF,logZMF] = UGM_Infer_MeanField(nodePot,edgePot,edgeStruct);
%     [MuEst,SigmaEst] = em_max(nodeBelMF,X(:));
%     MuEst
%     SigmaEst
%     diff = max(norm(MuEst-MuEst_o),norm(SigmaEst-SigmaEst_o));
%     i = i+1;
%   end
%   edgeStruct.maxIter = mf_iter_max;
%   solveMRF
% end


%% Run Chan Vese
I = X;
% Customerized Mask
m = zeros(size(I,1),size(I,2));
m(10:size(I,1)-10,10:size(I,2)-10) = 1;
mask = make_sdfunc(m);
prior_phis = mask; % just to get the size
seg = chanvese_exp(I,mask,600,0.1, [], 0, prior_phis,0,Epsilon,1,1); 

% % Built-in Mask
% seg = chanvese(I,'large',600,0.1); 
