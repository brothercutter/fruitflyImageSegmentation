addpath(genpath('../UGM/'))
addpath(genpath('../MRF'))

%% Algorithm constants
lambda = 1;
beta = 5;
Epsilon = 1;
superpixel = 1;

%% Load Data
% Load prior
load('../data/prior/shapePriorShift.mat');
% Load image
loadDataFruitfly;
% generateData;


%% Prior Pre-processing
Phi(Phi>0.5) = 1;
Phi(Phi<=0.5) =0;
prior_phis = Phi(:,:,[19:27]);
%prior_phis = Phi(:,:,[1,20]); 
no_ex = size(prior_phis,3);
% Compute kernel sigma
kernel_sigma = mean(get_Hdistance(circshift(prior_phis,[0 0 -1]),prior_phis,Epsilon));
%kernel_sigma = 1;

% Make priors and images match for CV
prior_size = size(Phi(:,:,1));
I = X;
I = mean(I,3);
I = (255 - I)/255;
X_CV = im2double(imresize(I,prior_size));
%[nRows,nCols] = prior_size;


%% Superpixels

% Create superpixels
I = im2double(X); 
createSuperPixels;
nSp = length(unique(Sp(:)));
[nRows, nCols] = size(I);
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

% Resize superpixels to ChanVese
Sp_r = imresize(Sp,prior_size,'nearest');
G = zeros(nSp,2);  % index matrix for centers of gravity in CV image size
G_vec = zeros(nSp,1);
no_sp = zeros(nSp,1);
for i = 1:nSp
    [row_i col_i] = find(Sp_r == i);
    ind{i} = find(Sp_r == i);
    no_sp(i) = length(ind{i});
    G(i,:) = ceil(mean([row_i col_i],1));
    G_vec(i) = ceil(mean(ind{i}));
end


% For speeding up:
% Calculate the diffmat between the first and the rest. 
% Then only need to compute diffmat betwen current and the first
%diff_phi_mat = Heaviside(repmat(prior_phis(:,:,1), [1 1 no_ex-1]), Epsilon) - Heaviside(sd_phis, Epsilon);


%% Level set function initialization
m = zeros(prior_size);
m(10:size(X_CV,1)-10,5:size(X_CV,2)-10) = 1;
% m0 = Phi(:,:,36);
% m(m0>0.5) = 1;
% m(m0<=0.5) = 0;
phi = m; %(:); %make_sdfunc(m);

%% Algorithm
% MRF Initialization
MuEst = [0.1,.9];
SigmaEst = [0.2,0.1];
%probA = ones(nNodes,nStates)*(1/nStates);
i = 0;
diff = 1;
iter_max = 20;
nodeBelMF_CV = zeros(prior_size);
 
for i = 1:5
    
  % Convert into MRF Superpixel stuff
  phi_vec = phi(:);
  expLambdaPhi = exp(lambda*phi_vec(G_vec));
  probLogistic = 1./(expLambdaPhi + 1);
  probA = [probLogistic, 1 - probLogistic]; 
  
%   
%   expLambdaPhi = exp(lambda*phi(:));
%   probLogistic = 1./(expLambdaPhi + 1);
%   probA = [probLogistic, 1 - probLogistic]; 
  
%   figure;
%   imagesc(reshape(probA(:,2)./probA(:,1),nRows,nCols));
%   title('prob(y_i|a_i)');
%   colorbar
%   
  %pause

  %nodePot_old = nodePot;
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
  imagesc(labelSuperpixelsOnImage(Sp,nodeBelMF(:,1)));
  colormap gray
  title('Mean Field Estimates of Marginals');
  fprintf('(paused)\n');

  % Convert into CV 
  for j = 1:nSp;      
       nodeBelMF_CV(ind{j},:) = repmat(nodeBelMF(1,:),[no_sp(j) 1]);
  end
  EQ_M = lambda*(probLogistic.*nodeBelMF_CV(:,2) - (1-probLogistic).*nodeBelMF_CV(:,1));
  EQ_M = reshape(EQ_M,size(X_CV,1),size(X_CV,2));
  m = chanvese_exp(X_CV,m,150,0.1,beta*EQ_M,1,prior_phis, kernel_sigma, Epsilon); 
  phi = make_sdfunc(m);   % phi is a matrix!
  %pause
end

     
