%%
function [IDXT,DSTT]=knnsearch(varargin)  
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% knnsearch: Linear k-nearest neighbor (KNN) search   在R中找Q的k邻域
% Q,R中每一行表示一个点的分量
% [IDX,DST] = knnsearch(Q,R,k) searches the reference data set R (L x d
% array representing L points in a d-dimensional space) to find 
% the k-nearest neighbors of each query point represented by each row of Q
% (m x d array). The results are stored in the (m x k) index array, IDX.
% The corresponding distances are stored in another array, DST.
%
% IDX = knnsearch(Q,R) takes the default value k=1.
% IDX = knnsearch(Q) or IDX = knnsearch(Q,[],k) does the search for R = Q.
%
% Modified from KNNSEARCH.m by Yi Cao at Cranfield University
% See also, kdtree, nnsearch, delaunary, dsearch, dsearchn
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Parse inputs
[Q,R,k,fident] = parse_inputs(varargin{:});
% Check outputs
error(nargoutchk(0,2,nargout));
%%
[m,d] = size(Q);
L=size(R,1);
IDXT = zeros(m,k);
DSTT = zeros(m,k);
%%
if k==1
    % Loop for each query point
    for i=1:m
        DSALL=zeros(L,1);
        for t=1:d
            DSALL=DSALL+(R(:,t)-Q(i,t)).^2;
        end
        if fident
            DSALL(i)=inf;
        end
        [DSTT(i),IDXT(i)]=min(DSALL);
    end
else
    for i=1:m
        DSALL=zeros(L,1);
        for t=1:d
            DSALL=DSALL+(R(:,t)-Q(i,t)).^2;
        end
        if fident
            DSALL(i)=inf;
        end
        [sorted_dst,sorted_idx]=sort(DSALL);
        IDXT(i,:)=sorted_idx(1:k);
        DSTT(i,:)=sorted_dst(1:k);
    end
end
%%
if nargout>1
    DSTT=sqrt(DSTT);
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [Q,R,k,fident] = parse_inputs(varargin)
%% Parse the inputs
error(nargchk(1,3,nargin));
%%
Q=varargin{1};
%%
if nargin<2
    fident = true;
    R=Q;
else
    fident = false;
    R=varargin{2};
end
%%
if isempty(R)
    fident = true;
    R=Q;
end
%%
if ~fident
    fident = isequal(Q,R);
end
%%
if nargin<3
    k=1;
else
    k=varargin{3};
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%