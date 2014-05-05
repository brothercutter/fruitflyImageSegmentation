
generateData;
probA = ones(1,nNodes)*(1/nStates);
edgeStruct = createImgStructure(X,nStates);
knownMuSigma = 0;
if knownMuSigma  
  [nodePot,edgePot] = producePotentials(X,edgeStruct,Mu,Sigma,probA);
  solveMRF;
else
  MuEst = [0.1,0.5];
  SigmaEst = [0.1,0.1];
  for i = 1:10
    [nodePot,edgePot] = producePotentials(X,edgeStruct,MuEst,SigmaEst,probA);
    [nodeBelMF,edgeBelMF,logZMF] = UGM_Infer_MeanField(nodePot,edgePot,edgeStruct);
    [MuEst,SigmaEst] = em_max(nodeBelMF,X(:));
    MuEst
    SigmaEst
  end
  solveMRF
end
