function varargout=evaluate_functions(xcell,varargin)
% INTERNAL FUNCTION: Evaluates functions in various formats
%
% ::
%
%   varargout=evaluate_functions(xcell,varargin)
%
% Args:
%
%    - **xcell** [cell|struct|fhandle]: function of interest
%
%      - input is a cell array: each element in a cell is assumed to be a
%        function handle
%      - input is a struct:
%
%        - if the field names are among: 'size','functions','map','partitions'
%          then we are dealing with derivatives
%        - if the field names are among: 'code','argins','argouts' then the
%          function could not be written as a function handle and has to be
%          evaluated using eval
%
%      - input is a function handle : direct evaluation
%
%    - **varargin** : input arguments to the function
%
% Returns:
%    :
%
%    - **varargout** : ouput arguments to the function
%


nout=nargout;

varargout=cell(1,nout);

if iscell(xcell)
    
    varargout{1}=main_engine();
    
    varargout{1}=sparse(cell2mat(varargout{1}));
    % varargout{1}=sparse(vec(cell2mat(varargout{1}.'))); % <-- 

elseif isstruct(xcell)
    
    derivative_fields={'size','derivatives'};% field 'map' removed to allow for vectorization
    
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
        
        vector_it=@(x)cell2mat(x(:).');
        
        for iout=1:nout
            
            xcell=tmp(iout).derivatives;
            
            if ~iscell(xcell)
                
                % xcell={xcell};                
                error('Unexpected format for derivatives');
       
            end
            
            locs=xcell(:,2);
            
            is_vectorized=size(locs{1},2)==1;
            
            if is_vectorized
                
                locs=locs{1};
                
            end
            
            xcell=xcell(:,1);
            
            mm=tmp(iout).size(1);
            
            nn=tmp(iout).size(2);
            
            if isempty(xcell)
                
                xout=sparse(mm,nn);
                
            else
                
                B=utils.code.set_vectorized_mapping(locs);
                
                vals=main_engine();
                
                nzmax=tmp(iout).nnz_derivs;
                
                ii=B(:,1);
                
                jj=B(:,2);
                
                v=vector_it(vals);
                
                v=v(B(:,3));
                
                xout=sparse(ii,jj,v,mm,nn,nzmax);
                
            end
            
            varargout{iout}=xout;
            
        end
        
    end

    function xout=main_engine()
        
        n=numel(xcell);
        
        xout=xcell;
        
        for item=n:-1:1
            
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
if ~isempty(code__)
    
    eval(code__)
    
end

% recover the outputs
%--------------------
varargout=argouts__(1:nout);

for ivar=1:nout
    
    varargout{ivar}=sparse(eval(varargout{ivar}));

end

end