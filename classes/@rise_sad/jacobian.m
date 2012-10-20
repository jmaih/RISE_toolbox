function [Jac,wrt,myrefs,i_index]=jacobian(objectives,varnames,wrt,myrefs,i_index)
% objectives is a function handle or an array of function handles
% wrt is a rise_sad variable or a cell array of such.
% output is returned in a rise_sad form in case the user wants to use it to
% recompute the hessian, which would be much easier to compute from the
% rise_sad form
if nargin<5
    i_index=[];
    if nargin<4
        myrefs=[];
        if nargin<3
            wrt=[];
            if nargin<2
                error([mfilename,':: the list of the arguments entering the function must be provided'])
            end
        end
    end
end

% probably possible to link different trees through some kind of hashing
% but I do not have that skill yet
%  example
%  func={'exp(a+b*log(c)+c*atan(a*b))','cos(abs(a)+atan(a*b))'};
%  [Jac,wrt,myrefs,i_index]=rise_sad.jacobian(func,{'a','c','b'})
%  
if isempty(i_index),i_index=0;end

if isempty(myrefs),myrefs={};end

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
%% locate the wrt variables in the args
wrt_locs=locate_variables(wrt,args_);
%% save the variable names and create the objects
for ii=1:numel(args)
    args{ii}=rise_sad(args{ii});
end
% variables going into the jacobian
wrt_args=args(wrt_locs);
%% recreate the functions
ncols=numel(wrt);
nrows=numel(objectives);
Jac=repmat({rise_sad(0)},nrows,ncols);
for irow=1:nrows
    [occur,myfunc]=find_occurrences(objectives{irow},args_);
    loc=ismember(wrt,args_(occur));
    % re-create the function
    argfun=cell2mat(strcat(args_(occur),','));
    myfunc=str2func(['@(',argfun(1:end-1),')',myfunc]);
    % create the master tree
    if any(occur)
        objectives{irow}=myfunc(args{occur});
    else
        objectives{irow}=rise_sad(0);
    end
    % differentiate the master tree
    if any(loc)
        Jac(irow,loc)=objectives{irow}.diff(wrt_args(loc));
    end
end
%% re_flag the master tree after the derivatives have been computed
for irow=1:nrows
    i_index=objectives{irow}.re_flag_tree(i_index);
end
%% collect the references
for irow=1:nrows
    myrefs=[myrefs;collect_references(objectives{irow})];
end
end

function [occur,objectives]=find_occurrences(objectives,vlist)
if isa(objectives,'function_handle')
    objectives=func2str(objectives);
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
    if isa(names{ii},'rise_sad')
        names{ii}=names{ii}.name;
    end
    if ~ischar(names{ii})
        error([mfilename,':: variable names must be char'])
    end
    bad=names{ii}(~isstrprop(names{ii},'alphanum'));
    if ~isempty(bad)
        error([mfilename,':: ',bad,' is(are) not valid character(s) is variable names '])
    end
end
for ii=1:numel(names)
    ni=sum(strcmp(names{ii},names));
    if ni>1
        error([mfilename,':: variable ',names{ii},' declared more than once'])
    end
end
end