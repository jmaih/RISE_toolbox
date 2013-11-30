function varargout=online_function_evaluator(F,varargin)

% evaluates 3 types of functions:
% A) function handles
% B) structures containing functions in cell arrays as well as other
% parameters determining the size output matrices and the location of the
% different elements computed in those matrices.
% C) string codes with inputs and outputs as follows:
% inputs:
% - F: structure with fields code, argins and argouts
% - varargin instances of elements whose names are given in argins
% the code may have references to nargin_ or nargout_ those variables will
% automatically be replaced by nargin and nargout in the formal function,
% while they will be evaluated in the code
% example of usage
% F=struct();
% F.code='a=1;b=a+p(5);c=diag([a,b]);Jac=c;';
% F.argins='p';
% F.argouts='Jac';
% test=online_function_evaluator(F,rand(10,1))
% see also: code2file

% code=struct('code',code,'argins',{{'y','x','ss','param','def'}},'argouts',{{'Jac'}})
Default=struct('code','char','argins','char or cellstr','argouts','char or cellstr');

if nargin==0
    varargout=Default;
    return
end

nargout_=nargout;
nargin_=nargin; %#ok<NASGU>
varargout=cell(1,nargout_);
if isempty(F)
    % don't do any thing
elseif isa(F,'function_handle')
    % function handles are ready to go
    %---------------------------------
    % this is just pure beauty
    [varargout{1:nargout_}]=(F(varargin{:}));
elseif isfield(F,'functions')
    is_mapped=isfield(F,'map');
    for iout=1:nargout_
        siz=F(iout).size;
        if isnan(siz)
            tmp=F(iout).functions{1}(varargin{:});
        else
            tmp=zeros(siz);%<--tmp=spalloc(siz(1),siz(2),F(iout).nnz_derivs);
            for irow=1:siz(1)
                thisfunc=F(iout).functions{irow};
                if ~isempty(thisfunc)
                    if is_mapped
                        tmp(irow,F(iout).map{irow})=thisfunc(varargin{:});
                    else
                        tmp(irow,:)=thisfunc(varargin{:});
                    end
                end
            end
        end
        varargout{iout}=sparse(tmp);
    end
else
    % here the code has to be evaluated
    %----------------------------------
    MajorFields=fieldnames(Default);
    for ifield=1:numel(MajorFields)
        thisfield=MajorFields{ifield};
        if isfield(F,thisfield)
            Default.(thisfield)=F.(thisfield);
        else
            Default.(thisfield)='';
        end
    end
    
    if isempty(Default.code)
        error('must have code')
    end
    if ~isempty(Default.argins) && ischar(Default.argins)
        Default.argins=cellstr(Default.argins);
    end
    if ~isempty(Default.argouts) && ischar(Default.argouts)
        Default.argouts=cellstr(Default.argouts);
    end
    
    n=length(varargin);
    if length(Default.argins)~=n
        error('number of inputs does not match the references in the code')
    end
    for ii=1:n
        param_i=Default.argins{ii};
        if ~ischar(param_i)
            error('names in argins should be char')
        end
        eval([param_i,'=varargin{ii};'])
    end
    
    eval(Default.code)
    
    varargout=Default.argouts(1:nargout_);
    for ivar=1:nargout_
        varargout{ivar}=sparse(eval(varargout{ivar}));
    end
end