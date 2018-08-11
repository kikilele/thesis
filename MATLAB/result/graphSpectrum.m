function s=graphSpectrum(adj) 
[~,D]=eig(laplacianMatrix(adj));
s=sort(diag(D)); % sort in decreasing order
s=s(2)/length(adj);