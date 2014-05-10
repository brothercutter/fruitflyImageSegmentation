% input Sp superpixel
addpath('../utilities/');
spCoord = zeros(nSp,2);

% plot the superpixel index on the segmented image
%{
figure; imagesc(I_sp);

for i = 1:nSp
  [y,x] = find(Sp == i);
  spCoord(i,:) = [mean(x),mean(y)];
end

text(spCoord(:,1), spCoord(:,2), num2strBatch(1:nSp)); 
hold on;
for e = 1:edgeStruct.nEdges
  n1 = edgeStruct.edgeEnds(e,1);
  n2 = edgeStruct.edgeEnds(e,2);

  plot(spCoord([n1,n2],1), spCoord([n1,n2],2));
  
end    

%}

se = strel('disk',1);
idxArray = cell(1,nSp);
for i = 1:nSp
  Y = zeros(nRows,nCols);
  Y(find(Sp == i)) = 1;
  dilatedY = imdilate(Y,se);
  idxArray{i} = find(dilatedY == 1);
  %figure;imagesc(dilatedY - Y); pause
end

adj = sparse(nSp,nSp);
for i = 1:(nSp-1)
  for j = (i+1):nSp
    idx1 = idxArray{i};
    idx2 = idxArray{j};
    if length(intersect(idx1,idx2)) > 20
      adj(i,j) = 1;
      adj(j,i) = 1;
    end
  end
end

edgeStruct = UGM_makeEdgeStruct(adj,nStates);