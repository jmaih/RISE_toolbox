classdef sad %< handle
    % to do: unknown functions
    % jacobian
    % hessian
    properties
        x
        dx
    end
    methods
        function obj=sad(y,dy)
            if nargin
                if isa(y,'sad')
                    obj=y;
                    return
                end
                if nargin<2
                    dy='0';
                end
                if ~ischar(y)||~ischar(dy)
                    error([mfilename,':: inputs must be char'])
                end
                obj.x=y;
                obj.dx=dy;
            end
        end
        function obj=plus(u,v)
            [u,du]=get_props(u);[v,dv]=get_props(v);
            xx=myplus(u,v);
            dxx=myplus(du,dv);
            obj=sad(xx,dxx);
        end
        function obj=uplus(u)
            [u,du]=get_props(u);
            xx=u;
            dxx=du;
            obj=sad(xx,dxx);
        end
        function obj=minus(u,v)
            [u,du]=get_props(u);[v,dv]=get_props(v);
            xx=myminus(u,v);
            dxx=myminus(du,dv);
            obj=sad(xx,dxx);
        end
        function obj=uminus(u)
            [u,du]=get_props(u);
            xx=myminus('0',u);
            dxx=myminus('0',du);
            obj=sad(xx,dxx);
        end
        function obj=times(u,v)
            obj=mtimes(u,v);
        end
        function obj=mtimes(u,v)
            [u,du]=get_props(u);[v,dv]=get_props(v);
            xx=mytimes(u,v);
            dxx=myplus(mytimes(du,v),mytimes(dv,u));
            obj=sad(xx,dxx);
        end
        function obj=power(u,v)
            obj=mpower(u,v);
        end
        function obj=mpower(u,v)
            [u,du]=get_props(u);[v,dv]=get_props(v);
            xx=mypower(u,v);
            dxx1=mytimes(mytimes(dv,['log(',u,')']),xx);
            dxx2=mytimes(mytimes(v,du),mypower(u,myminus(v,'1')));
            dxx=myplus(dxx1,dxx2);
            obj=sad(xx,dxx);
        end
        function obj=rdivide(u,v)
            obj=mrdivide(u,v);
        end
        function obj=mrdivide(u,v)
            [u,du]=get_props(u);[v,dv]=get_props(v);
            xx=mydivide(u,v);
            dxx=mytimes(du,v);
            dxx=myminus(dxx,mytimes(dv,u));
            dxx=mydivide(dxx,mypower(v,'2'));
            obj=sad(xx,dxx);
        end
        function obj=ldivide(u,v)
            obj=mldivide(u,v);
        end
        function obj=mldivide(u,v)
            obj=mrdivide(v,u);
        end
        function obj=exp(u)
            [u,du]=get_props(u);
            xx=['exp(',u,')'];
            dxx=mytimes(du,xx);
            obj=sad(xx,dxx);
        end
        function obj=log(u)
            [u,du]=get_props(u);
            xx=['log(',u,')'];
            dxx=mydivide(du,u);
            obj=sad(xx,dxx);
        end
        function obj=log10(u)
            obj=log(u)/log(10);
        end
        function obj=cos(u)
            [u,du]=get_props(u);
            xx=['cos(',u,')'];
            dxx=mytimes(myminus('0',du),['sin(',u,')']);
            obj=sad(xx,dxx);
        end
        function obj=acos(u)
            [u,du]=get_props(u);
            val=['acos(',u,')'];
            der=mypower(u,'2');
            der=myminus('1',der);
            der=mydivide(['-',parenthesize(du)],['sqrt(',der,')']);
            obj=sad(val,der);
        end
        function obj=cosh(u)
            [u,du]=get_props(u);
            val=['cosh(',u,')'];
            der=mytimes(du,['sinh(',u,')']);
            obj=sad(val,der);
        end
        function obj=sin(u)
            [u,du]=get_props(u);
            xx=['sin(',u,')'];
            dxx=mytimes(du,['cos(',u,')']);
            obj=sad(xx,dxx);
        end
        function obj=asin(u)
            [u,du]=get_props(u);
            val=['asin(',u,')'];
            der=mypower(u,'2');
            der=myminus('1',der);
            der=mydivide(du,['sqrt(',der,')']);
            obj=sad(val,der);
        end
        function obj=sinh(u)
            [u,du]=get_props(u);
            val=['sinh(',u,')'];
            der=mydivide(du,['cosh(',u,')']);
            obj=sad(val,der);
        end
        function obj=tan(u)
            [u,du]=get_props(u);
            val=['tan(',u,')'];
            der=mypower(['cos(',u,')'],'2');
            der=mydivide(du,der);
            obj=sad(val,der);
        end
        function obj=atan(u)
            [u,du]=get_props(u);
            val=['atan(',u,')'];
            der=mypower(u,'2');
            der=myplus('1',der);
            der=mydivide(du,['sqrt(',der,')']);
            obj=sad(val,der);
        end
        function obj=tanh(u)
            [u,du]=get_props(u);
            val=['tanh(',u,')'];
            der=mypower(['cosh(',u,')'],'2');
            der=mydivide(du,der);
            obj=sad(val,der);
        end
        function obj=min(u,v)
            if nargin~=2
                error([mfilename,':: number of arguments should be 2'])
            end
            [u,du]=get_props(u);[v,dv]=get_props(v);
            uLv=['(',u,'<',v,')'];
            val=['min(',u,',',v,')'];
            der=[uLv,'*',parenthesize(du),'+(1-',uLv,')*',parenthesize(dv)];
            obj = sad(val,der);
        end
        function obj=max(u,v)
            if nargin~=2
                error([mfilename,':: number of arguments should be 2'])
            end
            [u,du]=get_props(u);[v,dv]=get_props(v);
            uLv=['(',u,'<',v,')'];
            val=['max(',u,',',v,')'];
            der=[uLv,'*',parenthesize(dv),'+(1-',uLv,')*',parenthesize(du)];
            obj = sad(val,der);
        end
        function obj=sum(u,v)
            if nargin==1
                [u0,du0]=get_props(u(1));
                val=u0;
                der=du0;
                for ii=2:numel(u)
                    [u0,du0]=get_props(u(ii));
                    val=myplus(val,u0);
                    der=myplus(der,du0);
                end
                obj=sad(val,der);
            else
                obj=plus(u,v);
            end
        end
        function obj=normpdf(u,mu,sig)
            if nargin<3
                sig='1';
                if nargin<2
                    mu='0';
                end
            end
            [u,du]=get_props(u);
            mu=get_props(mu);
            sig=get_props(sig);
            val=['normpdf(',u,',',mu,',',sig,')'];
            der0=myminus(u,mu);
            der1=mypower(sig,'2');
            der=mydivide(['-(',der0,')'],der1);
            der=mytimes(der,du);
            der=mytimes(der,val);
            obj = sad(val,der);
        end
        function obj=normcdf(u,mu,sig)
            if nargin<3
                sig='1';
                if nargin<2
                    mu='0';
                end
            end
            [u,du]=get_props(u);
            mu=get_props(mu);
            sig=get_props(sig);
            val=['normcdf(',u,',',mu,',',sig,')'];
            der=mytimes(du,['normpdf(',u,',',mu,',',sig,')']);
            obj = sad(val,der);
        end
        function obj=abs(u)
            [u,du]=get_props(u);
            val=['abs(',u,')'];
            der=mytimes(['sign(',u,')'],du);
            obj=sad(val,der);
        end
        function obj=isreal(u)
            u=get_props(u);
            val=['real(',u,')'];
            der='0';
            obj=sad(val,der);
        end
        function obj=sqrt(u)
            [u,du]=get_props(u);
            xx=['sqrt(',u,')'];
            dxx=mytimes('0.5',mydivide(du,xx));
            obj=sad(xx,dxx);
        end
        function obj=norm(u)
            obj = sqrt(sum(u.^2));
        end
        
        function d = char(u,expand)
            if nargin<2
                expand=false;
            end
            d = u.dx;
            if expand
                d=sad.replace_keys(d);
            end
        end
    end
    methods(Static)
        varargout=jacobian(varargin)
        varargout=hessian(varargin)
        function [derivatives,auxiliary]=trim(derivatives)
            mapObj_sad=evalin('base','mapObj_sad');
            allKeys = keys(mapObj_sad);
            allValues = values(mapObj_sad);
            auxiliary='';
            % for some versions of matlab, one has to double the output of
            % the container, which returns elements of type int64 on my
            % machine: an undesirable feature.
            for icount=double(mapObj_sad.Count):-1:1 % icount=mapObj_sad.Count:-1:1
                % if a particular value occurs only once, then replace it
                % and remove the key
                vv=['T_',int2str(icount)];
                pat=['(?<![\w])',vv,'(?![\w])'];
                locs=regexp(derivatives,pat);
                if iscell(locs)
                    locs=cell2mat(locs(:)');
                end
                howmany=numel(locs);
                loc=strcmp(vv,allValues);
                add_it=true;
                if howmany==1
                    % normally I should do this for each cell,
                    % carefully reconstructing the string,... but this
                    % is not without cost.
                    derivatives=regexprep(derivatives,pat,['(',allKeys{loc},')']);
                    % do not add it if it does not already appear in the
                    % auxilary
                    if isempty(auxiliary)||isempty(regexp(auxiliary,pat,'once'))
                        add_it=false;
                    end
                end
                if add_it
                    auxiliary=[[vv,'=',allKeys{loc},';'],auxiliary];  %#ok<AGROW>
                end
                allKeys(loc)=[];
                allValues(loc)=[];
            end
            if ~isempty(auxiliary)
                auxiliary=strcat(regexp(auxiliary,';','split'),';');
                auxiliary=auxiliary(1:end-1);
            end
        end
        function string=replace_keys(string)
            % get the list of keys in the string
            mapObj_sad=evalin('base','mapObj_sad');
            allKeys = keys(mapObj_sad);
            allValues = values(mapObj_sad);
            pat='(?<![\w])T_[\d]+(?![\w])';
            string0=string;
            [start,finish]=regexp(string0,pat,'start','end');
            % reconstructing the string
            while ~isempty(start)
                string=string0(1:start(1)-1);
                for iloc=1:numel(start)
                    val=string0(start(iloc):finish(iloc));
                    val_loc= strcmp(val,allValues);
                    key=allKeys{val_loc};
                    parenth=(~isempty(string) && any(string(end)=='-^/*\'))||...
                        (finish(iloc)<length(string0) && any(string0(finish(iloc)+1)=='*/^\'));
                    if parenth
                        key=['(',key,')']; %#ok<AGROW>
                    end
                    suffix='';
                    if iloc<numel(start)
                        suffix=string0(finish(iloc)+1:start(iloc+1)-1);
                    end
                    string=[string,key,suffix]; %#ok<AGROW>
                end
                string=[string,string0(finish(end)+1:end)]; %#ok<AGROW>
                string0=string;
                [start,finish]=regexp(string0,pat,'start','end');
            end
        end
        function destroy()
            evalin('base','clear(''mapObj_sad'')');
        end
    end
end

function [u,du]=get_props(x)
du='0';
switch class(x)
    case 'sad'
        u=x.x;
        du=x.dx;
    case 'char'
        u=x;
    case 'double'
        u=num2str(x,10);
    otherwise
        error([mfilename,':: unsupported class ',class(x)])
end
[u,du]=archive(u,du);
end

function c=myplus(a,b)
% this operation is commutative and so sort alphabetically before computing
% this will help in the archivation process to make sure a+b=b+a
[a,b]=commute(a,b);
if strcmp(a,'0')
    if strcmp(b,'0')
        c='0';
    else
        c=tryevaluate(b);
    end
else
    if strcmp(b,'0')
        c=tryevaluate(a);
    else
        c=tryevaluate([a,'+',b]);
    end
end
end

function c=myminus(a,b)
if strcmp(a,'0')
    if strcmp(b,'0')
        c='0';
    else
        c=tryevaluate(['-',parenthesize(b,'+-')]);
    end
else
    if strcmp(b,'0')
        c=a;
    else
        c=tryevaluate([a,'-',parenthesize(b,'+-')]);
    end
end
end

function c=mytimes(a,b)
% this operation is commutative and so sort alphabetically before
% computing. this will help in the archivation process to make sure a*b=b*a
[a,b]=commute(a,b);

if strcmp(a,'0')||strcmp(b,'0')
    c='0';
elseif strcmp(a,'1')
    c=b;
elseif strcmp(b,'1')
    c=a;
else
    c=tryevaluate([parenthesize(a,'+-'),'*',parenthesize(b,'+-')]);
end
end

function c=mypower(a,b)
if strcmp(b,'0')||strcmp(a,'1')
    c='1';
else
    c=tryevaluate([parenthesize(a,'+-/*^'),'^',parenthesize(b,'+-/*^')]);
end
end

function c=mydivide(a,b)
if strcmp(a,'0')
    c='0';
else
    c=tryevaluate([parenthesize(a,'+-'),'/',parenthesize(b,'+-/*^')]);
end
end

function [a,b]=commute(a,b)
% sort alphabetically
if a(1)>b(1)
    atmp=a;
    btmp=b;
    a=btmp;
    b=atmp;
end
end

function x=parenthesize(x,forbid)
if nargin<2
    forbid='+-*/^';
end
forbid=strrep(forbid,'-','\-'); % must escape the minus sign
flag=~isempty(regexp(x,['[',forbid,']'],'start'));
if flag
    x=['(',x,')'];
end
end

function flag=is_atom(string)
flag=isempty(regexp(string,'[/*\-+^]','start','once')); %<---flag=isempty(regexp(string,'[/*\-+^]','start'));
end

function a=tryevaluate(a)
% checks whether a string can be evaluated
flag=~any(isstrprop(a,'alpha'));
if flag
    cntrl=a;
    cntrl(isstrprop(cntrl,'digit'))=[];
    flag=~isempty(cntrl) && ~isequal(cntrl,'.');
    if flag
        flag=false;
        for ii=1:length(cntrl)
            if any(cntrl(ii)=='+-*^/')
                flag=true;
                break
            end
        end
        if flag
            a=num2str(eval(a),10);
        end
    end
end
end

function varargout=archive(varargin)
if evalin('base','exist(''mapObj_sad'',''var'')')
    mapObj_sad=evalin('base','mapObj_sad');
else
    mapObj_sad=containers.Map;
end
varargout=varargin;
for item=1:length(varargin)
    if ~is_atom(varargin{item})
        if isKey(mapObj_sad,varargin{item})
            varargout{item}=mapObj_sad(varargin{item});
        else
            Count=mapObj_sad.Count+1;
            varargout{item}=['T_',int2str(Count)];
            mapObj_sad(varargin{item})=varargout{item};
        end
    end
end
assignin('base','mapObj_sad',mapObj_sad)
end

% function varargout=archive(varargin)
% persistent mapObj_sad
% if isempty(mapObj_sad)
%     mapObj_sad=cell(0,3);
% end
% varargout=varargin;
% for item=1:length(varargin)
%     if ~is_atom(varargin{item})
%         key=get_key(varargin{item});
%         if ~isempty(key) % <---if isKey(mapObj_sad,varargin{item})
%             varargout{item}=key;
%         else
%             if isempty(mapObj_sad)
%                 Count=1;
%             else
%                 Count=mapObj_sad{end,3}+1; % <--- Count=mapObj_sad.Count+1;
%             end
%             varargout{item}=['T_',int2str(Count)];
%             mapObj_sad=[mapObj_sad;varargin(item),varargout(item),{Count}]; %#ok<AGROW> % <-- mapObj_sad(varargin{item})=varargout{item};
%         end
%     end
% end
% % capture x and the derivatives right here. Check whether
% % they look like anything that has been computed earlier in
% % which case they should assume the auxiliary values... The
% % criterion should be that they have at least one operator
% % in order to be eligible for the container
% assignin('base','mapObj_sad',mapObj_sad)
%     function key=get_key(item)
%         key=[];
%         key_loc=find(strcmp(item,mapObj_sad(:,1)));
%         if ~isempty(key_loc)
%             key=mapObj_sad{key_loc,2};
%         end
%     end
% end