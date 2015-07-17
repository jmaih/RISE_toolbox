function xout=listing(varargin)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

% LISTING parses the arguments of a parameter or a variable
myscalar=@(x)isa(x,'double') && isscalar(x) && floor(x)==ceil(x);
default_xout={
    'name','',@(x)isvarname(x)% all
    'current_name','',@(x)isvarname(x)% all
    'tex_name','',@(x)ischar(x)% all
    'is_in_use',false,@(x)islogical(x)% all exogenous and parameters
    'governing_chain','',@(x)myscalar(x) && x>0% all parameters
    'is_switching',false,@(x)islogical(x)% all parameters
    'is_measurement_error',false,@(x)islogical(x)% all parameters
    'is_trans_prob',false,@(x)islogical(x)% all parameters
    'is_log_var',false,@(x)islogical(x)% all endogenous
    'is_endogenous',false,@(x)islogical(x)% all observables
    'max_lead',0,@(x)myscalar(x) && x>=0% all variables and parameters
    'max_lag',0,@(x)myscalar(x) && x<=0% all variables and parameters
    'state_id',nan,@(x)myscalar(x) && x>=0% all observables
	'is_auxiliary',false,@(x)islogical(x)% all endogenous
    };
    ff=default_xout(:,1).';
    nff=numel(ff);
    xout0=cell(1,nff);
[xout0{1:nff}]=utils.miscellaneous.parse_arguments(default_xout,varargin{:});

xout=struct();
for ifield=1:nff
    xout.(ff{ifield})=xout0{ifield};
end

if nargin
    if isempty(xout.name)
        error('must have a name')
    end
    if isempty(xout.current_name)
        xout.current_name=xout.name;
    end
else
    xout=xout(1:0);
end

end
