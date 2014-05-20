knownMuSigma = 1;
generateData;
probA = ones(1,nNodes)*(1/nStates);
edgeStruct = createImgStructure(X,nStates);
if knownMuSigma  
  [nodePot,edgePot] = producePotentials(X,edgeStruct,Mu,Sigma,probA);
  solveMRF;
else
  

