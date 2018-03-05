function d=dirichlet_shortcuts(a,location,lpdfn,rndfn)
% dirichlet_shortcuts -- memoizer for some dirichlet routines
%
% Syntax
% -------
% ::
%
%   d=dirichlet_shortcuts()
%
%   d=dirichlet_shortcuts(a,location)
%
%   d=dirichlet_shortcuts(a,location,lpdfn,rndfn)
%
% Inputs
% -------
%
% - **a** [k x 1 vector]: of hyperparameters for the distribution
%
% - **location** [k-1 x 1]: location of the parameters to estimate. Note
% that this vector has k-1 elements instead of k
%
% - **lpdfn** [function_handle]: function handle for the log pdf of the
% dirichlet distribution
%
% - **rndfn** [function_handle]: function handle for random draws of the
% dirichlet distribution
%
% Outputs
% --------
%
% - **d** [struct]: with the following elements
%   - **lpdfn** [function_handle]: takes as input an k-1 x 1 vector but
%   returns the log density for an k x 1 vector
%   - **rndfn** [function_handle]: returns and k-1 x n matrix of random
%   draws of the dirichlet distribution
%   - **location** [vector]: location of the parameters of interest in the
%   vector of estimated parameters
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

if nargin==0
    
    d=struct('lpdfn',{},'rndfn',{},'location',{});
    
else
            
        if nargin<4
            
            rndfn=[];
            
        end
                
    if isempty(rndfn)
        
        [lpdfn,~,~,rndfn]=distributions.dirichlet();
        
    end
    
    d=struct('lpdfn',@logdensity,'rndfn',@dirich_draw,'location',location);
    
end

    function lpdf=logdensity(x)
        % add the last guy since RISE does not include all the elements in
        % the vector.
        theta=[x(:);1-sum(x)];
        
        lpdf=lpdfn(theta,a);
        
    end

    function d=dirich_draw(n)
        
        if nargin==0||isempty(n)
            
            n=1;
            
        end
        
        d=rndfn(a,[],n);
        % keep only the relevant since RISE does not include all the
        % elements in the vector
        d=d(1:end-1);
        
    end

end
