% input Sp superpixel
addpath('../utilities/');
spCoord = zeros(nSp,2);

% plot the superpixel index on the segmented image
%{
for i = 1:nSp
  [y,x] = find(Sp == i);
  spCoord(i,:) = [mean(x),mean(y)];
end

text(spCoord(:,1), spCoord(:,2), num2strBatch(1:nSp)); 
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
    if length(intersect(idx1,idx2)) > 0
      adj(i,j) = 1;
      adj(j,i) = 1;
    end
  end
end

edgeStruct = UGM_makeEdgeStruct(adj,nStates);