function fv=fiedlerVector(adj)  
[V,D]=eig(laplacianMatrix(adj));
[~,Y]=sort(diag(D));
fv=V(:,Y(2));