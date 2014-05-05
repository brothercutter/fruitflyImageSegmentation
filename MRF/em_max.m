function [Mu, Sigma] = em_max(dist, img)

% dist: probability distribution matrix, size nxl, 
% img: observed image intensities % vector

n = length(img);  % number of pixels
l = size(dist,2); % number of labels

Mu = zeros(1,l);
Sigma = zeros(1,l);
for k = 1:l
  w = dist(:,k);
  [mu,sigma2] = weightMeanVar(img,w);
  Mu(k) = mu;
  Sigma(k) = sqrt(sigma2);
end


% Calculate number of pixels per label
%no_pix = dist*ones(n,l);  % size lx1

% Calculate mean
%mu = dist*img./no_pix; % size lx1

% Calculate variance
%var_vec = (repmat(img,[l 1]) - repmat(mu, [1 n])).^2;  % size lxn
%sigma = dist'*var_vec'./no_pix;

