addpath(genpath('../UGM/'))
addpath(genpath('../MRF'))

%% Algorithm constants
lambda = 1;
beta = 1;
Epsilon = 1;
superpixel = 1;

% Load prior
load('../data/prior/shapePriorShift.mat');
%% Prior Pre-processing
Phi = 1 - Phi;
Phi(Phi>0.5) = 1;
Phi(Phi<=0.5) =0;
prior_phis = Phi(:,:,19:24);
no_ex = size(prior_phis,3);
% Compute kernel sigma
kernel_sigma = mean(get_Hdistance(circshift(prior_phis,[0 0 -1]),prior_phis,Epsilon));

% loadData
%I = imread('../data/embryo.png');
loadDataFruitfly;
I = X;
I = im2double(I);
[nRows,nCols,dump] = size(I);

% Make priors and images match for CV
prior_size = size(Phi(:,:,1));
I_CV = imresize(mean(I,3),prior_size);
I_CV = 1 - I_CV;
I_CV = I_CV/max(I_CV(:));


rng(215)
%createSuperPixels;
I = 1-mean(I,3);
I = I/max(I(:));
nSp = length(unique(Sp(:)));

% Resize superpixels to ChanVese
Sp_r = imresize(Sp,prior_size,'nearest');
nSp = length(unique(Sp_r));
idx0 = unique(Sp_r);
for i = 1:nSp
    Sp_r(find(Sp_r == idx0(i))) = i;
end  
Sp = imresize(Sp_r,size(I(:,:,1)),'nearest');
nStates = 2;
createNeighbor;


G = zeros(nSp,2);  % index matrix for centers of gravity in CV image size
no_sp = zeros(nSp,1);
for i = 1:nSp
    [row_i col_i] = find(Sp_r == i);
    ind{i} = find(Sp_r == i);
    no_sp(i) = length(ind{i});
    G(i,:) = ceil(mean([row_i col_i],1));    
end
figure; imagesc(Sp_r);
text(G(:,2), G(:,1), num2strBatch(1:nSp)); 


% transform the pixel intensity to superpixel intensity
X = zeros(nSp,1);
for i = 1:nSp
  X(i) = mean(I(find(Sp == i)));
end
%X(9) = 0;
figure;
imagesc(labelSuperpixelsOnImage(Sp,X));
colorbar;



%% Level set function initialization
m = zeros(prior_size);
m(10:size(I_CV,1)-10,5:size(I_CV,2)-10) = 1;
% m0 = Phi(:,:,36);
% m(m0>0.5) = 1;
% m(m0<=0.5) = 0;
phi = m; %(:); %make_sdfunc(m);

%% Algorithm
% MRF Initialization
MuEst = [0.1,.9];
SigmaEst = [0.1,0.2];
%probA = ones(nNodes,nStates)*(1/nStates);

i = 0;
diff = 1;
iter_max = 20;
nodeBelMF_CV = zeros(prior_size(1)*prior_size(2),nStates);
probLogistic_CV =zeros(prior_size(1)*prior_size(2),1);
for i = 1:5
    
  % Convert into MRF Superpixel stuff
  phi_vec = phi(:);
  idxTemp = G(:,1)+(G(:,2)-1)*prior_size(1);
  expLambdaPhi = exp(lambda*phi(idxTemp));
  
  probLogistic = 1./(expLambdaPhi + 1);
  
  
  probA = [probLogistic, 1 - probLogistic]; 
  
   figure;
   imagesc(labelSuperpixelsOnImage(Sp_r,probA(:,2)));
   title('prob(y_i=1|a_i)');
   colorbar


   [MuEst,SigmaEst] = normalMixtureEM(MuEst,SigmaEst,X,edgeStruct,probA,5);
   [nodePot,edgePot] = producePotentials(X,edgeStruct,MuEst,SigmaEst,probA);
   nodePot_norm = sum(nodePot,2);
  
%   figure;
%   imagesc(reshape(nodePot(:,2)./nodePot_norm,nRows,nCols));
%   title('Node Potential Difference');  
%   colorbar
%   pause
  
  [nodeBelMF,edgeBelMF,logZMF] = UGM_Infer_MeanField(nodePot,edgePot,edgeStruct);  
  
  figure;
  imagesc(labelSuperpixelsOnImage(Sp,nodeBelMF(:,2)));
  colorbar;
  title('Mean Field Estimates of Marginals (y=2)');
  fprintf('(paused)\n');
  pause

  % Convert into CV 
  for j = 1:nSp;      
       nodeBelMF_CV(ind{j},:) = repmat(nodeBelMF(j,:),[no_sp(j) 1]);
       probLogistic_CV(ind{j}) = probLogistic(j); 
     end
  EQ_M = lambda*(probLogistic_CV.*nodeBelMF_CV(:,2) - (1-probLogistic_CV).*nodeBelMF_CV(:,1));
  EQ_M = reshape(EQ_M,size(I_CV,1),size(I_CV,2));
  %m = chanvese_exp(I_CV,m,150,0.1,beta*EQ_M,1,prior_phis, kernel_sigma, Epsilon); 
  m = chanvese_exp(I_CV,m,150,0.1,nodeBelMF_CV,1,prior_phis, kernel_sigma, Epsilon,beta,lambda); 
  
  phi = make_sdfunc(m);   % phi is a matrix!
  pause
end

     
