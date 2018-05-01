function [newstruct,retcode]=evaluate_transition_matrices(xstruct,varargin)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:


if isa(xstruct,'function_handle')
    
    [newstruct,retcode]=xstruct(varargin{:});
    
elseif isstruct(xstruct)
    
    newstruct=struct();
    
    chain_names=fieldnames(xstruct);
    
    is_loose_commit=any(strcmp(chain_names,parser.loose_commit()));
    
    if is_loose_commit
        
        Qinit=1;
        
    end
    
    Q=1;
    
    for iname=1:numel(chain_names)
        
        chain=chain_names{iname};
        
        if ~isa(xstruct.(chain),'function_handle')
            
            error('all elements in xstruct should be function handles')
            
        end
        
        newstruct.(chain)=xstruct.(chain)(varargin{:});
        
        Q=kron(Q,newstruct.(chain));
        
        if is_loose_commit && ~strcmp(chain,parser.loose_commit())
            
            Qinit=kron(Qinit,newstruct.(chain));
            
        end
        
    end
    
    [newstruct.Q,retcode]=utils.code.validate_transition_matrix(Q);
    
    if ~retcode
        
        if ~is_loose_commit
            
            Qinit=Q;
            
        end
        
        [newstruct.Qinit,retcode]=utils.code.validate_transition_matrix(Qinit);
        
    end
    
else
    
    error('first input argument must be a structure')
    
end

end