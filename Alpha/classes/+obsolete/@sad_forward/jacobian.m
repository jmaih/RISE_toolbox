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
prototype=sad_forward('xxx','0');

%% save the variable names and create the objects
args_= args;
for w=1:numel(args)
    prototype.x=args{w};
    args_{w}=prototype;
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
    [occur,myfunc]=find_occurrences(objectives{irow},args);
    
    % re-create the function
    var_occur=args(occur);
    argfun=cell2mat(strcat(var_occur,','));
    myfunc=str2func(['@(',argfun(1:end-1),')',myfunc]);
    
    % wrt variables occurring in var_occur SILENTLY
    locs=locate_variables(wrt,var_occur,true);
    wrt_occur=wrt(~isnan(locs));
    locs_=locs(~isnan(locs));
    nwrt=numel(wrt_occur);
    
    % now inflate the derivatives of the guys occurring
    var_occur_=args_(occur);
    dx={'0'};
    dx=dx(1,ones(nwrt,1));
    for ii=1:nwrt
        dxi=dx;
        dxi{ii}='1';
        var_occur_{locs_(ii)}.dx=dxi;
    end
    
    % locate the specific wrt
    orig_locs=locate_variables(wrt_occur,wrt);
    if ~isempty(orig_locs)
        % differentiate
        tmp=myfunc(var_occur_{:});
        Jac(irow,orig_locs)=char(tmp);
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