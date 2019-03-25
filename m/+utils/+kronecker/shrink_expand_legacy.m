%  INTERNAL FUNCTION: computes shrinking and expansion objects for the manipulation of symmetric tensors
% 
%  ::
% 
%    [keep,expand,C,UC,B]=shrink_expand(n,k)
%    [keep,expand,C,UC,B]=shrink_expand(n,k,strategy)
%    [keep,expand,C,UC,B]=shrink_expand(n,k,strategy,debug)
% 
%  Args:
% 
%     - **n** [scalar] : number of variables in the tensor
%     - **k** [scalar] : order of the tensor
%     - **strategy** [1|2|{3}] : alternative computation strategies
% 
%       - 1: uses bsxfun
%       - 2: splanar inspired
%       - 3: applies ismember
% 
%     - **debug** [true|{false}] : checks the results
% 
%  Returns:
%     :
% 
%     - **keep** [logical] : n^k x 1 vector true for the columns to be kept
%     - **expand** [vector] : n^k x 1 vector replicating the compressed
%       columns to form the grand tensor
%     - **C** [matrix] : sparse compression matrix of size n^k x g, where
%       g=nchoosek(n+k-1,k) is the number of unique elements in the tensor
%       (matrix version of **keep**)
%     - **UC** [matrix] : sparse expansion matrix of size g x n^k, where
%       g=nchoosek(n+k-1,k) is the number of unique elements in the tensor
%       (matrix version of **expand**)
%     - **B** [matrix] : g x k matrix of combinations without repetitions. Each
%       row in increasing order.
% 
%  Note:
% 
%     useful for shrinking tensors of the form fvvv...v as used in higher-order
%     differentiation
% 
%