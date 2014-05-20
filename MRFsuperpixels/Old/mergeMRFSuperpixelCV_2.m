close all

addpath(genpath('../UGM/'))
addpath(genpath('../MRF'))

%% Algorithm constants
lambda = 1;
beta = 1.5;
Epsilon = 10;
superpixel = 1;

% Load prior
load('../data/prior/shapePriorShift.mat');

%% Prior Pre-processing
Phi = 1 - Phi;
Phi(Phi>0.5) = 1;
Phi(Phi<=0.5) =0;
prior_phis = Phi(:,:,1:27);
no_ex = size(prior_phis,3);
% Compute kernel sigma
kernel_sigma = mean(get_Hdistance(circshift(prior_phis,[0 0 -1]),prior_phis,Epsilon));
%kernel_sigma=1e-2;
no_ex = size(prior_phis,3);

% Create sd_phis
for k = 1:no_ex
    sd_phis(:,:,k) = make_sdfunc(prior_phis(:,:,k));
end


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

%% Superpixels
rng(215)
createSuperPixels;
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
%figure; imagesc(Sp_r);
%text(G(:,2), G(:,1), num2strBatch(1:nSp));
idxTemp = G(:,1)+(G(:,2)-1)*prior_size(1);

% transform the pixel intensity to superpixel intensity
X = zeros(nSp,1);
for i = 1:nSp
    X(i) = mean(I(find(Sp == i)));
end

figure;
imagesc(labelSuperpixelsOnImage(Sp,X));
colorbar;



%% Level set function initialization
m = zeros(prior_size);
m(10:size(I_CV,1)-10,5:size(I_CV,2)-10) = 1;
%  m0 = Phi(:,:,23);
%  m(m0>0.5) = 1;
%  m(m0<=0.5) = 0;
phi = make_sdfunc(m); %m; %(:); %

%% Algorithm
% MRF Initialization
%MuEst = [0.4,.5]; %0.41, 0.75
%SigmaEst = [0.13,0.119];

MuEst = [0.1,.9];
SigmaEst = [0.1,0.1];

i = 0;
diff = 1;
iter_max = 20;
nodeBelMF_CV = zeros(prior_size(1)*prior_size(2),nStates);
probLogistic_CV =zeros(prior_size(1)*prior_size(2),1);

size_ICV = size(I_CV);
phi_vec = phi(:);
expLambdaPhi = exp(lambda*phi(idxTemp));
probLogistic = 1./(expLambdaPhi + 1);
probA = ones(nSp,nStates)*(1/nStates);

for i = 1:10
    
    
    %probA = ones(nSp,nStates)*(1/nStates);
    
    figure;
    imagesc(labelSuperpixelsOnImage(Sp_r,probA(:,2)));
    title('prob(y_i=1|a_i)');
    colorbar
    
    
    [MuEst,SigmaEst] = normalMixtureEM(MuEst,SigmaEst,X,edgeStruct,probA,20);
    [nodePot,edgePot] = producePotentials(X,edgeStruct,MuEst,SigmaEst,probA);
    nodePot_norm = sum(nodePot,2);
    [nodeBelMF,edgeBelMF,logZMF] = UGM_Infer_MeanField(nodePot,edgePot,edgeStruct);
    
    figure;
    imagesc(labelSuperpixelsOnImage(Sp,nodeBelMF(:,2))); colorbar;
    title('Mean Field Estimates of Marginals (y=1)');
    %pause
    
    % Convert into CV
    for j = 1:nSp;
        nodeBelMF_CV(ind{j},:) = repmat(nodeBelMF(j,:),[no_sp(j) 1]);
        %probLogistic_CV(ind{j}) = probLogistic(j);
    end
%     EQ_M = lambda*(probLogistic_CV.*nodeBelMF_CV(:,2) - (1-probLogistic_CV).*nodeBelMF_CV(:,1));
%     EQ_M = reshape(EQ_M,size_ICV(1),size_ICV(2));
    %m = chanvese_exp(I_CV,m,150,0.1,beta*EQ_M,1,prior_phis, kernel_sigma, Epsilon);
    
    % Take the most likely binary mask of MRF as initialization
    
%     mask = getmask_frompot(nodeBelMF_CV,size_ICV(1),size_ICV(2));
%     mask = make_sdfunc(mask);
    mask = phi;
    m = chanvese_exp(I_CV,mask,150,0.1,nodeBelMF_CV,1,sd_phis, kernel_sigma, Epsilon,beta,lambda);
    
    phi = make_sdfunc(m);   % phi is a matrix!
    % Convert into MRF Superpixel stuff
    phi_vec = phi(:);
    expLambdaPhi = exp(lambda*phi(idxTemp));
    probLogistic = 1./(expLambdaPhi + 1);
    probA = [probLogistic, 1 - probLogistic];
    %pause
end


