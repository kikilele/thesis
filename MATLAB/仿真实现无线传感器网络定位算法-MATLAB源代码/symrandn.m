% Generate a Normal random matrix of size k, but with zero diagonal
% and ensure the matrix is symmetric (noise(i,j) == noise(j,i)).
function [noise] = symrandn(k)

temp = tril(randn(k),-1);
noise = temp + temp';