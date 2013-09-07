function [jaco_tree,code,code_expanded]=jacobian(objectives,varlist,wrt,generic,start)
% objectives is a function handle or an array of function handles
% wrt is a sadiff variable or a cell array of such.
if nargin<5
    start=[];
    if nargin<4
        generic='';
    end
end
if isempty(start)
    start=1;
end
if isempty(generic)
    generic='Jac';
end

is_setup=isa(objectives,'sadiff') && isa(wrt,'sadiff');
wrt_index={};
if is_setup
    mytree=objectives;
else
    [mytree,wrt_index,wrt]=sadiff.setup(objectives,varlist,wrt);
end

jaco_tree=diff(mytree,wrt);

if nargout>1
    [c,mycall]=print(jaco_tree,wrt_index);
    indx_prefix=mycall.prefix_list{end};
    feed=strcat(mycall.fid(:,1),'=',mycall.fid(:,2),';');
    code=replace_strings(feed,c,generic,start,indx_prefix);
    if nargout>2
        discard_useless_parentheses=true;
        code_expanded=char(jaco_tree,[],[],discard_useless_parentheses);
        % code_expanded
        if ischar(code_expanded)
            code_expanded={code_expanded};
        end
        nn=start-1+(1:numel(code_expanded));
        nn=num2str(nn');
        code_expanded=strcat(generic,'(',nn,',:)=',code_expanded(:),';');
    end
end


function fid=replace_strings(fid,c,generic,start,indx_prefix)
if nargin<5
    indx_prefix=[];
    if nargin<4
        start=1;
        if nargin<3
            generic='Jac';
        end
    end
end
indx_name=':';
for idef=1:numel(c)
    if ~isempty(indx_prefix)
    indx_name=create_handle(idef,indx_prefix);
    end
    if any(c{idef}=='_') % cheap way of recognizing valid definitions
        fid=strrep(fid,c{idef},...
            [generic,'(',sprintf('%0.10g',start+idef-1),...
            ',',indx_name,')']);
    end
end

%
% [~,dxhandle,~,~,mycall] = differentiate(mytree,wrt);
%
% eqtn_nbr=numel(mytree);
% wrt_nbr=numel(wrt);
% JJ2=strcat(mycall.fid(:,1),'=',mycall.fid(:,2),';');
% for idef=1:numel(dxhandle)
%     JJ2=strrep(JJ2,dxhandle{idef},['Jac(',sprintf('%0.10g',idef),',:)']);
% end
% JJ2=[['Jac=zeros(',int2str(eqtn_nbr),',',int2str(wrt_nbr),');']
%     ['bigi_=speye(',int2str(wrt_nbr),');']
%     JJ2];
