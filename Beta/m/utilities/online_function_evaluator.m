function varargout=online_function_evaluator(F,varargin)

% evaluates string code as a function with inputs and outputs
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
if isa(F,'function_handle')
%     % get the number of output arguments of the function
%     nout=nargout(F);
%     varargout=cell(1,nout);
%     % this is just pure beauty
%     [varargout{1:nout}]=(F(varargin{:}));
    varargout=cell(1,nargout_);
    [varargout{1:nargout_}]=(F(varargin{:}));
else
    MajorFields=fieldnames(Default);
    for ifield=1:numel(MajorFields)
        thisfield=MajorFields{ifield};
        if isfield(F,thisfield)
            Default.(thisfield)=F.(thisfield);
            F=rmfield(F,thisfield);
        else
            Default.(thisfield)='';
        end
    end
    remaining_fields=fieldnames(F);
    if ~isempty(remaining_fields)
        disp(remaining_fields)
        error('these fields are not recognized')
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
        varargout{ivar}=eval(varargout{ivar});
    end
end