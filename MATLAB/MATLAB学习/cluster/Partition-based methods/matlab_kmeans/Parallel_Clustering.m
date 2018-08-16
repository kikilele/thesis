% Randomly generate a large data set from a Gaussian mixture model.
clear;
Mu = bsxfun(@times,ones(20,30),(1:20)'); % Gaussian mixture mean
rn30 = randn(30,30);
Sigma = rn30'*rn30; % Symmetric and positive-definite covariance 对称和正定协方差
Mdl = gmdistribution(Mu,Sigma);

rng(1); % For reproducibility
X = random(Mdl,10000);

% Mdl is a 30-dimensional gmdistribution model with 20 components. 
% X is a 10000-by-30 matrix of data generated from Mdl.


pool = parpool;                      % Invokes workers
stream = RandStream('mlfg6331_64');  % Random number stream
options = statset('UseParallel',1,'UseSubstreams',1,...
    'Streams',stream);

tic; % Start stopwatch timer
[idx,C,sumd,D] = kmeans(X,20,'Options',options,'MaxIter',10000,...
    'Display','final','Replicates',10);
toc % Terminate stopwatch timer






