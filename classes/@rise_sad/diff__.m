function [derivs,wrt] = diff(func,wrt)
% rise_sad/DIFF overloads diff with a rise_sad object argument
% DIFF(S,'v')
% DIFF(S,'v1','v2',...,'vn')
% DIFF(S,'v1',m1,'v2',m2,...,'vn',mn)
% 0- S is a string or an anonymous function on which one can
% apply func2str
% 1- detect the variables list
% 2- check that the wrt variables are present. if they are not
% then they automatically get derivative 0
create_func=false;
if ischar(func)
    fstring=func;
    create_func=true;
elseif isa(func,'function_handle')
    fstring=func2str(func);
else
    error([mfilename,':: first argument must be a string or a function handle'])
end

% get the symbolic expression and the list of the variables present in the
% function
validNames=valid_varnames();
[fstring,varlist]=analytical_symbolic_form(fstring,validNames,'symbolic');
args=strcat(varlist,',');
args=[args{:}];
func_declare=['@(',args(1:end-1),')'];
if create_func
    func=[func_declare,fstring];
    func=str2func(func);
end
if nargin<2
    wrt=varlist;
end
if ischar(wrt)
    wrt=cellstr(wrt);
end
nargs=length(wrt);
% check that the names are valid
wrtnames=cell(1,0);
for ivar=1:nargs
    tmp=wrt{ivar};
    if ischar(tmp)
        wrtnames=[wrtnames,{tmp}]; %#ok<*AGROW>
    elseif iscellstr(tmp)
        wrtnames=[wrtnames,tmp(:)'];
    else
        error([mfilename,':: arguments should be cellstr or char'])
    end
end
% check that the names are valid
valid_varnames(wrtnames);

% initialize the variables and their derivatives
nvars=numel(varlist);
x=repmat({rise_sad('xxx','0')},1,nvars);
for ivar=1:nvars
    x{ivar}.x=varlist(ivar);
end
nder=numel(wrt);
derivs=repmat({'0'},1,nder);
silence=true;
for ider=1:nder
    locs=locate_variables(wrt{ider},varlist,silence);
    if all(~isnan(locs)) % all variables must be present else the derivative is zero
        func_ider=func;
        for iloc=1:numel(locs)
            z=x;
            z{locs(iloc)}.dx={'1'};
            derivs{ider}=char(func_ider(z{:}));
            original=derivs{ider};
            if isempty(derivs{ider})
                keyboard
            end
            derivs{ider}=rise_sad.optimize(derivs{ider});
            if isempty(derivs{ider})
                keyboard
            end
            if iloc<numel(locs)
                func_ider=str2func([func_declare,derivs{ider}]);
            end
        end
        derivs{ider}=analytical_symbolic_form(derivs{ider},validNames,'analytical');
    end
end
end
