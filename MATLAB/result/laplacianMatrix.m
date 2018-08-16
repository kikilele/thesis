function L=laplacianMatrix(adj) 
L=diag(sum(adj))-adj;