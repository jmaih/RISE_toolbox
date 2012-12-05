function [Jac,auxiliary]=jacobian(objectives,varnames,wrt)
% objectives is a function handle or an array of function handles
% wrt is a sad_forward variable or a cell array of such.
% output is returned in a sad_forward form in case the user wants to use it to
% recompute the hessian, which would be much easier to compute from the
% sad_forward form
    if nargin<3
        wrt=[];
        if nargin<2
            error([mfilename,':: the list of the arguments entering the function must be provided'])
        end
    end

    
if ischar(objectives)
    objectives=cellstr(objectives);
end
if ~iscell(objectives)
    objectives={objectives};
end

if ischar(varnames)
    varnames=cellstr(varnames);
end
args=get_names(varnames);
if isempty(wrt)
    wrt=args;
end
wrt=get_names(wrt);

% locate the variables in main list
bad=~ismember(wrt,args);
bad=wrt(bad);
if ~isempty(bad)
    disp(bad)
    error([mfilename,':: the variables above not present in the list of variables'])
end

args_=args;

%% save the variable names and create the objects
args=repmat({sad_forward('xxx')},numel(args),1);
for ii=1:numel(args_)
    args{ii}.x=args_{ii};
end
%% recreate the functions
ncols=numel(wrt);
nrows=numel(objectives);
Jac=repmat({'0'},nrows,ncols);
for irow=1:nrows
    if isa(objectives{irow},'sad_forward')
        % expand it in full
        objectives{irow}=char(objectives{irow},1);
    end
    [occur,myfunc]=find_occurrences(objectives{irow},args_);
    
    % re-create the function
    var_occur=args_(occur);
    argfun=cell2mat(strcat(var_occur,','));
    myfunc=str2func(['@(',argfun(1:end-1),')',myfunc]);
    
    % differentiate
    for icol=1:numel(wrt)
        args_x=args(occur);
        target= strcmp(wrt{icol},var_occur);
        if any(target)
            args_x{target}.dx='1';
            Jac{irow,icol}=myfunc(args_x{:});
            Jac{irow,icol}=char(Jac{irow,icol});
        end
    end
end

auxiliary={};
end

function [occur,objectives]=find_occurrences(objectives,vlist)
if isa(objectives,'function_handle')
    objectives=func2str(objectives);
    if strcmp(objectives(1),'@')
        first_close=find(objectives==')',1,'first');
        objectives=objectives(first_close+1:end);
    end
end
if ~ischar(objectives)
    error([mfilename,':: first input must be a string or a function handle'])
end
if ischar(vlist)
    vlist=cellstr(vlist);
end
vlist=vlist(:)';
args=cell2mat(strcat(vlist,'|'));
varlist = regexp(objectives,['(?<![\w])(',args(1:end-1),')(?![\w])'],'match');
occur=ismember(vlist,varlist);
end

function names=get_names(array)
if ~iscell(array)
    array={array};
end
names=array;
for ii=1:numel(names)
    if isa(names{ii},'sad_forward')
        names{ii}=names{ii}.name;
    end
    if ~ischar(names{ii})
        error([mfilename,':: variable names must be char'])
    end
    bad=regexp(names{ii},'[\W]','match');%<-- names{ii}(~isstrprop(names{ii},'alphanum'));
    if ~isempty(bad)
        bad=cell2mat(bad);
        error([mfilename,':: ''',bad,''' is(are) not valid character(s) is variable names '])
    end
end
for ii=1:numel(names)
    ni=sum(strcmp(names{ii},names));
    if ni>1
        error([mfilename,':: variable ',names{ii},' declared more than once'])
    end
end
end