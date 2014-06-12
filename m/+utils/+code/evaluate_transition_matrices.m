function [newstruct,retcode]=evaluate_transition_matrices(xstruct,varargin)

newstruct=struct();
if ~isstruct(xstruct)
    error('first input argument must be a structure')
end
chain_names=fieldnames(xstruct);
if ~all(isfield(xstruct,chain_names))
    error('chain_names should be fields of xstruct')
end

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

[newstruct.Q,retcode]=validate_transition_matrix(Q);

if ~retcode
    if ~is_loose_commit
        Qinit=Q;
    end
    [newstruct.Qinit,retcode]=validate_transition_matrix(Qinit);
end

    function [Qx,rcode]=validate_transition_matrix(Qx)
        total=sum(Qx,2);
        if any(isnan(Qx(:))) || ...
                any(Qx(:)<0) || ...
                any(Qx(:)>1)|| ...
                ~max(abs(total-1))<1e-9
            rcode=3;
        else
            rcode=0;
            Qx=bsxfun(@rdivide,Qx,total);
        end
    end
end