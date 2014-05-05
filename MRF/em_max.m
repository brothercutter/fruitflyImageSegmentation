function [mu, sigma] = em_max(dist, img)

% dist: probability distribution matrix, size nxl, 
% img: observed image intensities % vector

n = length(img);  % number of pixels
l = size(dist,1); % number of labels

% Calculate number of pixels per label
no_pix = dist'*ones(n,l);  % size lx1

% Calculate mean
mu = dist'*img./no_pix; % size lx1

% Calculate variance
var_vec = (repmat(img,[l 1]) - repmat(mu, [1 n])).^2;  % size lxn
sigma = dist'*var_vec'./no_pix;

