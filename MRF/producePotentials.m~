function [nodePot,edgePot] =  producePotentials(X,edgeStruct,Mu,Sigma,probA)
% X: the image;
% Mu: estimate for the means;
% Sigma: estimate for the standard deviations;
% probA: the estimated prob(y_i|a_i) vector
% create note and edge potentials
% probA = ones(1,nNodes)*(1/nStates);

Xvec = X(:);
Xstd = UGM_standardizeCols(Xvec,1);
nNodes = length(Xvec);
nStates = size(Mu,2);

nodePot = zeros(nNodes,nStates);

for i = 1:nNodes
  for l = 1:nStates
    nodePot(i,l) = (1/(sqrt(6.28)*Sigma(l)))*exp(-(Xvec(i)-Mu(l))^2/(2*Sigma(l)^2))*probA(i);
  end
end

edgePot = zeros(nStates,nStates,edgeStruct.nEdges);
for e = 1:edgeStruct.nEdges
  n1 = edgeStruct.edgeEnds(e,1);
  n2 = edgeStruct.edgeEnds(e,2);

  pot_same = exp(1.8 + .3*1/(1+abs(Xstd(n1)-Xstd(n2))));
  edgePot(:,:,e) = [pot_same 1;1 pot_same];
end