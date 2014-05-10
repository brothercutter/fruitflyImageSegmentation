function boundary_im = findboundary(X)

se = strel('disk',1);
dilatedIm = imdilate(X,se);
boundary_im = dilatedIm - X;
