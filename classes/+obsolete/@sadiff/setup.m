function [mytree,wrt_index,wrt,mycall]=setup(objectives,varlist,wrt)
% objectives is a function handle or an array of function handles
% wrt is a sadiff variable or a cell array of such.

% % [objectives_done,objectives]=check_is_done(objectives);

if ischar(varlist),varlist=cellstr(varlist);end

[wrt_done,wrt]=check_is_done(wrt);

varlist=get_names(varlist);

%% prototypical sadiff object
prototype=sadiff.prototype();

%% make the wrt variables into objects?
if ~wrt_done
    if isempty(wrt)
        wrt=varlist;
    end
    wrt=get_names(wrt);
    for irt=1:numel(wrt)
        func=wrt{irt};
        wrt{irt}=prototype;
        wrt{irt}.func=func;
        wrt{irt}.order=irt;
    end
    % put at a (row) vector
    wrt=[wrt{:}];
end

wrt_names={wrt.func};
%% save the variable names and create the objects
varlist_= varlist;

%% dealing with the objectives
if isa(objectives,'sadiff')
    objectives=num2cell(objectives);
end
objectives=objectives(:);
nobj=numel(objectives);
mytree=sadiff.empty(0);
wrt_index=cell(nobj,1);
for iobj=1:nobj
    [occur,thisFunc]=find_occurrences(objectives{iobj},varlist);
    var_occur=varlist(occur);
    wrt_index{iobj}=find(ismember(wrt_names,var_occur));
    if isempty(wrt_index{iobj})
        mytree(iobj,1)=sadiff(0);
    else
        if isa(thisFunc,'sadiff')
            mytree(iobj,1)=thisFunc;
        else
            if ~isa(varlist_{1},'sadiff')
                nvar=numel(varlist);
                for w=1:nvar
                    prototype.func=varlist{w};
                    varlist_{w}=prototype;
                end
            end
            % re-create the function
            argfun=cell2mat(strcat(var_occur,','));
            thisFunc=str2func(['@(',argfun(1:end-1),')',thisFunc]);
            % tree building
            obj_args=varlist_(occur);
            mytree(iobj,1)=thisFunc(obj_args{:});
        end
    end
end

%% create the global element
mycall=sadiff.metastruct();

end

function [objectives_done,objectives]=check_is_done(objectives)
objectives_done=false;
if ~isempty(objectives)
    switch class(objectives)
        case 'char'
            objectives=cellstr(objectives);
        case {'cellstr','cell'}
        case 'sadiff'
            objectives_done=true;
        otherwise
            error(['class ',class(objectives),' not allowed'])
    end
    if ~objectives_done && isa(objectives,'sadiff')
        objectives=[objectives{:}];
        objectives_done=true;
    end
end

end

function [occur,objectives]=find_occurrences(objectives,vlist)
if isa(objectives,'sadiff')
    varlist=load_varlist(objectives);
else
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
    varlist=cell2mat(strcat(vlist,'|'));
    varlist = regexp(objectives,['(?<![\w])(',varlist(1:end-1),')(?![\w])'],'match');
end
occur=ismember(vlist,varlist);

end

function names=get_names(array)
if ~iscell(array)
    array={array};
end
names=array;
for ii=1:numel(names)
    if isa(names{ii},'sadiff')
        names{ii}=names{ii}.func;
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
