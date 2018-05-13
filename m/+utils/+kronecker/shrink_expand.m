function [keep,expand,C,UC]=shrink_expand(n,k,varargin)
% shrink_expand computes shrinking and expansion objects for the manipulation of symmetric tensors
%
% ::
%
%
%   [keep,expand,C,UC,B]=shrink_expand(n,k)
%   [keep,expand,C,UC,B]=shrink_expand(n,k,strategy)
%   [keep,expand,C,UC,B]=shrink_expand(n,k,strategy)
%
% Args:
%
%    - **n** [scalar] : number of variables in the tensor
%    - **k** [scalar] : order of the tensor
%    - **strategy** ['up'|{'down'}] : alternative compression strategies
%      - 'up': i1<=i2<=i3....<=in
%      - 'down': i1>=i2>=i3....>=in
%
% Returns:
%    :
%
%    - **keep** [cell array of k logical vectors] : n^k x 1 vector true for
%    the columns to be kept  
%    - **expand** [cell array of k vectors] : n^k x 1 vector replicating the
%    compressed columns to form the grand tensor
%    - **C** [cell array of k matrices] : sparse compression matrix of size
%    n^k x g, where g=nchoosek(n+k-1,k) is the number of unique elements in
%    the tensor (matrix version of **keep**)
%    - **UC** [cell array of k matrices] : sparse expansion matrix of size
%    g x n^k, where g=nchoosek(n+k-1,k) is the number of unique elements in
%    the tensor (matrix version of **expand**)
%
% Note:
%
%    Essentially the results are given for orders from 1 to k
%
%    useful for shrinking tensors of the form fvvv...v as used in higher-order
%    differentiation
%
% Example:
%
%    See also:

nout=nargout;

[expand,keep]=utils.kronecker.tensor_locator(n,k,varargin{:});

if nout>2
    
    C=cell(1,k);
    
    UC=cell(1,k);
    
end

for ii=1:k
    
    expand{ii}=expand{ii}(:).';
    
    keep{ii}=sparse(keep{ii}(:).');
    
    matrix_form()
    
end

    function matrix_form()
        
        if nout>2
            
            nbig=numel(keep{ii});
            
            Ci=speye(nbig);
            
            C{ii}=Ci(:,keep{ii});
            
            if nout>3
                
                nkeep=sum(keep{ii});
                
                UCi=speye(nkeep);
                
                UC{ii}=UCi(:,expand{ii});
                
            end
            
        end
        
    end

end
