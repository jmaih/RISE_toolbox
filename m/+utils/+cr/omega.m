function oo = omega(nz,index)
% omega -- creator of sparse matrices for mutivariate chain rules up to
% fifth order
%
% ::
%
%
%   oo = omega(nz,index)
%
% Args:
%
%    - **nz** [integer]: number of variables in the differentiation
%
%    - **index** [integer]: code for the requested omega. Must be in [1,9]
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

% Reference: Oren Levintal(2014): "Fifth Order Perturbation Solution to
% DSGE Models", 2014. 

j =  2; k =  3; l =  4; m =  5; n =  6;

switch index
    
    case 1
        
        oo=isum([k,l,j],[j,l,k],[j,k,l]);
    
    case 2
        
        oo=isum([l,m,k,j],[k,m,l,j],[j,m,l,k],[k,l,m,j],...
            [j,l,m,k],[j,k,m,l]);
    
    case 3
        
        oo=isum([k,l,m,j],[j,l,m,k],[j,k,m,l],[j,k,l,m]);
    
    case 4
        
        oo=isum([k,l,j,m],[j,l,k,m],[j,k,l,m]);
    
    case 5
        
        oo=isum([m,n,l,k,j],[l,n,m,k,j],[k,n,m,l,j],...
            [j,n,m,l,k],[l,m,n,k,j],[k,m,n,l,j],...
            [j,m,n,l,k],[k,l,n,m,j],[j,l,n,m,k],...
            [j,k,n,m,l]);
    
    case 6
        
        oo=isum([l,m,n,k,j],[k,m,n,l,j],[j,m,n,l,k],...
            [k,l,n,m,j],[j,l,n,m,k],[j,k,n,m,l],...
            [k,l,m,n,j],[j,l,m,n,k],[j,k,m,n,l],...
            [j,k,l,n,m]);
    
    case 7
        
        oo=isum([l,m,k,n,j],[k,m,l,n,j],[j,m,l,n,k],...
            [k,l,m,n,j],[j,l,m,n,k],[j,k,m,n,l],...
            [l,m,j,n,k],[k,m,j,n,l],[j,m,k,n,l],...
            [k,l,j,n,m],[j,l,k,n,m],[j,k,l,n,m],...
            [k,l,j,m,n],[j,l,k,m,n],[j,k,l,m,n]);
    
    case 8
        
        oo=isum([k,l,m,n,j],[j,l,m,n,k],[j,k,m,n,l],...
            [j,k,l,n,m],[j,k,l,m,n]);
    
    case 9
        
        oo=isum([k,l,m,j,n],[j,l,m,k,n],[j,k,m,l,n],...
            [j,k,l,m,n],[k,l,n,j,m],[j,l,n,k,m],...
            [j,k,n,l,m],[j,m,n,k,l],[k,m,n,j,l],[l,m,n,j,k]);
    
    otherwise
        
        error('index must be in [1,9]')

end

    function ss=isum(varargin)
        
        order=numel(varargin{1});
        
        nzi=nz^order;
        
        M=reshape(1:nzi,[1,nz*ones(1,order)]);
        
        Iz=speye(nzi);
        
        ss=sparse(nzi,nzi);
        
        for ii=1:length(varargin)
            
            ss=ss+Iz(:,ipermute(M,[1,varargin{ii}]));
        
        end
        
    end

end

