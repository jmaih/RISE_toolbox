function varargout=evaluate_functions(xcell,varargin)

nout=nargout;
varargout=cell(1,nout);
if iscell(xcell)
    varargout{1}=main_engine();
    varargout{1}=sparse(cell2mat(varargout{1}));
elseif isstruct(xcell)
    derivative_fields={'size','functions','map','partitions'};%,'maxcols','nnz_derivs'
    eval_fields={'code','argins','argouts'};
    if all(isfield(xcell,derivative_fields))
        [varargout{1:nout}]=derivative_engine();
    elseif all(isfield(xcell,eval_fields))
        [varargout{1:nout}]=eval_engine(xcell,nout,varargin{:});
    else
        error('unknown type of structure')
    end
elseif isa(xcell,'function_handle')
    [varargout{1:nout}]=xcell(varargin{:});
else
    error('first input must be a cell, a structure or a function handle')
end

    function varargout=derivative_engine()
        tmp=xcell;
        varargout=cell(1,nout);
        for iout=1:nout
            xcell=tmp(iout).functions;
            xout=zeros(tmp.size);
            if ~isempty(xcell)
                vals=main_engine();
                for irows=1:size(xcell,1)
                    if ~isempty(vals{irows})
                        xout(irows,tmp(iout).map{irows})=vals{irows};
                    end
                end
            end
            varargout{iout}=sparse(xout(:,tmp(iout).partitions));
        end
    end

    function xout=main_engine()
        n=numel(xcell);
        xout=xcell;
        for item=1:n
            if ~isempty(xcell{item})
                xout{item}=xcell{item}(varargin{:});
            end
        end
    end
end

function varargout=eval_engine(xcell,nout,varargin)
argins__=xcell.argins;
argouts__=xcell.argouts;
code__=xcell.code;
if numel(argouts__)<nout
    error('too many output arguments')
end
if numel(argins__)~=length(varargin)
    error('wrong number of input arguments')
end
% evaluate input arguments
%-------------------------
for iarg=1:numel(argins__)
    param_i=argins__{iarg};
    if ~ischar(param_i)
        error('names in argins should be char')
    end
    eval([param_i,'=varargin{iarg};'])
end
% evaluate the code
%------------------
eval(code__)
% recover the outputs
%--------------------
varargout=argouts__(1:nout);
for ivar=1:nout
    varargout{ivar}=sparse(eval(varargout{ivar}));
end
end