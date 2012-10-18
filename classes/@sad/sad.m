% record_book={iter,var,expression}
% record_book=cell(0,3);
% derivatives go into the record book
% all constants have empty book
classdef sad
    properties
        x
        dx
    end
    methods
        function obj=sad(u,du)
            if nargin
                obj.x=u;
                obj.dx=du;
                if exist('sad_records','var')
                    sad.update_book(obj)
                else
                    assignin('base','sad_records',cell(0,2))
                end
            end
        end
        function obj=mtimes(u,v)
            [u,du]=get_props(u);[v,dv]=get_props(v);
            val=mytimes(u,v);
            der=myplus(...
                mytimes(du,v),...
                mytimes(dv,u)...
                );
            obj=sad.update_book(val,der);
        end
        function obj=times(u,v)
            obj=mtimes(u,v);
        end
        function obj=plus(u,v)
            [u,du]=get_props(u);[v,dv]=get_props(v);
            val=myplus(u,v);
            der=myplus(du,dv);
            obj=sad.update_book(val,der);
        end
        function obj=uplus(u)
            [u,du]=get_props(u);
            % no need to update the book
            obj=sad(u,du);
        end
        function obj=minus(u,v)
            [u,du]=get_props(u);[v,dv]=get_props(v);
            val=myminus(u,v);
            der=myminus(du,dv);
            obj=sad.update_book(val,der);
        end
        function obj=uminus(u)
            [u,du]=get_props(u);
            % no need to update the book
            obj=sad(['-(',u,')'],['-(',du,')']);
        end
        function obj=mrdivide(u,v)
            [u,du]=get_props(u);[v,dv]=get_props(v);
            val=mydivide(u,v);
            if strcmp(dv,'0')
                der=mydivide(du,v);
            else
                der=mytimes(du,v);
                der=myminus(der,mytimes(dv,u));
                der=mydivide(der,mypower(v,'2'));
            end
            obj=sad.update_book(val,der);
        end
        function obj=rdivide(u,v)
            obj=mrdivide(u,v);
        end
        function obj=mldivide(u,v)
            obj=mrdivide(v,u);
        end
        function obj=ldivide(u,v)
            obj=mldivide(u,v);
        end
        function obj=mpower(u,v)
            [u,du]=get_props(u);[v,dv]=get_props(v);
            val=mypower(u,v);
            if strcmp(dv,'0')
                der=mytimes(mytimes(v,du),mypower(u,myminus(v,'1')));
            else
                der1=mytimes(dv,['log(',u,')']);
                der2=mytimes(mydivide(du,u),v);
                der=mytimes(myplus(der1,der2),val);
            end
            obj=sad.update_book(val,der);
        end
        function obj=power(u,v)
            obj=mpower(u,v);
        end
        function obj=sqrt(u)
            obj=u^0.5;
        end
        function obj=log(u)
            [u,du]=get_props(u);
            val=['log(',u,')'];
            der=mydivide(du,u);
            obj=sad.update_book(val,der);
        end
        function obj=log10(u)
            obj=rdivide(log(u),'log(10)');
        end
        function obj=exp(u)
            [u,du]=get_props(u);
            val=['exp(',u,')'];
            der=mytimes(du,val);
            obj=sad.update_book(val,der);
        end
        function obj=cos(u)
            [u,du]=get_props(u);
            val=['cos(',u,')'];
            der=mytimes(['-sin(',u,')'],du);
            obj=sad.update_book(val,der);
        end
        function obj=acos(u)
            [u,du]=get_props(u);
            val=['acos(',u,')'];
            der=mypower(u,'2');
            der=myminus('1',der);
            der=mydivide(['-',parenthesize(du)],['sqrt(',der,')']);
            obj=sad.update_book(val,der);
        end
        function obj=cosh(u)
            [u,du]=get_props(u);
            val=['cosh(',u,')'];
            der=mytimes(du,['sinh(',u,')']);
            obj=sad.update_book(val,der);
        end
        function obj=sin(u)
            [u,du]=get_props(u);
            val=['sin(',u,')'];
            der=mytimes(['cos(',u,')'],du);
            obj=sad.update_book(val,der);
        end
        function obj=asin(u)
            [u,du]=get_props(u);
            val=['asin(',u,')'];
            der=mypower(u,'2');
            der=myminus('1',der);
            der=mydivide(du,['sqrt(',der,')']);
            obj=sad.update_book(val,der);
        end
        function obj=sinh(u)
            [u,du]=get_props(u);
            val=['sinh(',u,')'];
            der=mydivide(du,['cosh(',u,')']);
            obj=sad.update_book(val,der);
        end
        function obj=tan(u)
            [u,du]=get_props(u);
            val=['tan(',u,')'];
            der=mypower(['cos(',u,')'],'2');
            der=mydivide(du,der);
            obj=sad.update_book(val,der);
        end
        function obj=atan(u)
            [u,du]=get_props(u);
            val=['atan(',u,')'];
            der=mypower(u,'2');
            der=myplus('1',der);
            der=mydivide(du,['sqrt(',der,')']);
            obj=sad.update_book(val,der);
        end
        function obj=tanh(u)
            [u,du]=get_props(u);
            val=['tanh(',u,')'];
            der=mypower(['cosh(',u,')'],'2');
            der=mydivide(du,der);
            obj=sad.update_book(val,der);
        end
        function obj=abs(u)
            [u,du]=get_props(u);
            val=['abs(',u,')'];
            der=mytimes(['sign(',u,')'],du);
            obj=sad.update_book(val,der);
        end
        function obj=real(u)
            [u,du]=get_props(u);
            val=['real(',u,')'];
            der=['real(',du,')'];
            obj=sad.update_book(val,der);
        end
        function obj=norm(u)
            obj = sqrt(sum(u.^2));
        end
        function obj=sum(u,v)
            % rise_sad/SUM overloads sum with a rise_sad object argument
            if nargin==1
                [u0,du0]=get_props(u(1));
                val=u0;
                der=du0;
                for ii=2:numel(u)
                    [u0,du0]=get_props(u(ii));
                    val=myplus(val,u0);
                    der=myplus(der,du0);
                end
                obj=sad.update_book(val,der);
            else
                obj=plus(u,v);
            end
        end
        function obj=normpdf(u,mu,sig)
            if nargin<3
                sig=1;
                if nargin<2
                    mu=0;
                end
            end
            [u,du]=get_props(u);
            [mu,sig]=get_char_form(mu,sig);
            val=['normpdf(',u,',',mu,',',sig,')'];
            der0=myminus(u,mu);
            der1=mypower(sig,'2');
            der=mydivide(['-(',der0,')'],der1);
            der=mytimes(der,du);
            der=mytimes(der,val);
            obj = sad.update_book(val,der);
        end
        function obj=normcdf(u,mu,sig)
            if nargin<3
                sig=1;
                if nargin<2
                    mu=0;
                end
            end
            [u,du]=get_props(u);
            [mu,sig]=get_char_form(mu,sig);
            val=['normcdf(',u,',',mu,',',sig,')'];
            der=mytimes(du,['normpdf(',u,',',mu,',',sig,')']);
            obj = sad.update_book(val,der);
        end
        function obj=conj(u)
            [u,du]=get_props(u);
            val=['conj(',u,')'];
            der=['conj(',du,')'];
            obj = sad.update_book(val,der);
        end
        function c = isreal(u)
            % rise_sad/ISREAL overloads isreal with a rise_sad object argument
            u=get_props(u);
            c = ['isreal(',u,')'];
        end
        function obj = min(u,v)
            % sad/MIN overloads min with a rise_sad object argument
            % but will work with 2 arguments.
            if nargin~=2
                error([mfilename,':: number of arguments should be 2'])
            end
            [u,du]=get_props(u);[v,dv]=get_props(v);
            uLv=['(',u,'<',v,')'];
            val=['min(',u,',',v,')'];
            der=[uLv,'*',parenthesize(du),'+(1-',uLv,')*',parenthesize(dv)];
            obj = sad.update_book(val,der);
        end
        function obj = max(u,v)
            % sad/MAX overloads max with a rise_sad object argument
            % but will work with 2 arguments.
            if nargin~=2
                error([mfilename,':: number of arguments should be 2'])
            end
            [u,du]=get_props(u);[v,dv]=get_props(v);
            uLv=['(',u,'<',v,')'];
            val=['max(',u,',',v,')'];
            der=[uLv,'*',parenthesize(dv),'+(1-',uLv,')*',parenthesize(du)];
            obj = sad.update_book(val,der);
        end
        function d = char(u)
            d = u.dx;
        end
    end
    methods(Static)
        function obj=update_book(val,der)
            % Now we can also simplify equalities
            % check that the expression was not computed before
            sad_records=evalin('base','sad_records');
            olditer=size(sad_records,1);
            push_val=assign_element(val);
            push_der=assign_element(der);
            obj=sad(push_val,push_der);
            assignin('base','sad_records',sad_records)
            function push_element=assign_element(element)
                push_element=element;
                is_definition=any(strcmp(element,sad_records(:,1)));
                if ~strcmp(element,'0') && ~strcmp(element,'1') && ~is_definition
                    % is is possible to re-arrange the expression wrt the
                    % operators so as to facilitate the search?
                    % check whether the whole expression has been computed
                    % before.
                    % if not, check the bits of the expression that have
                    % been computed before
                    dejavu=find(strcmp(element,sad_records(:,2)));
                    if ~isempty(dejavu)
                        push_element=sad_records{dejavu,1};
                    else
                        iter=olditer+1;
                        push_element=['T_',int2str(iter)];
                        sad_records=[sad_records
                            {push_element,element}];
                        olditer=iter;
                    end
                end
            end
        end
        varargout=diff(varargin)
        varargout=hessian(varargin)
        varargout=jacobian(varargin)
    end
end

function varargout=get_char_form(varargin)
n=length(varargin);
for ivar=1:n
    vv=varargin{ivar};
    switch class(vv)
        case 'char'
            varargout{ivar}=vv;
        case 'double'
            varargout{ivar}=num2str(vv,10);
        case 'rise_sad'
            varargout{ivar}=vv.x;
        otherwise
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
end

function c=myplus(a,b)
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
        c=tryevaluate(['-',parenthesize(b)]);
    end
else
    if strcmp(b,'0')
        c=a;
    else
        c=tryevaluate([a,'-',parenthesize(b)]);
    end
end
end

function c=mytimes(a,b)
if strcmp(a,'0')||strcmp(b,'0')
    c='0';
elseif strcmp(a,'1')
    c=b;
elseif strcmp(b,'1')
    c=a;
else
    c=tryevaluate([parenthesize(a),'*',parenthesize(b)]);
end
end

function c=mypower(a,b)
if strcmp(b,'0')
    c='1';
else
    c=tryevaluate([parenthesize(a),'^',parenthesize(b)]);
end
end

function c=mydivide(a,b)
if strcmp(a,'0')
    c='0';
else
    c=tryevaluate([parenthesize(a),'/',parenthesize(b)]);
end
end

function x=parenthesize(x)
flag=false;
for ii=1:length(x)
    if any(x(ii)=='+-*/^')
        flag=true;
        break
    end
end
if flag
    x=['(',x,')'];
end
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
