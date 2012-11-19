function c=rise_cat(a,b,func)

switch func
    case {'plus','minus'}
        if strcmp(a,'0')
            if strcmp(b,'0')
                c='0';
            else
                if strcmp(func,'minus')
                    b=['-',parenthesize(b,'+-')];
                end
                c=b;
            end
        else
            if strcmp(b,'0')
                c=a;
            else
                if strcmp(func,'minus')
                    if ~rise_isa(b,'atom') && ~rise_isa(b,'function')
                        b=parenthesize(b,'+-');
                    end
                    c=[a,'-',b];
                else
                    [a,b]=commute(a,b);
                    c=[a,'+',b];
                end
            end
        end
    case 'times'
        if strcmp(a,'0')||strcmp(b,'0')
            c='0';
        elseif strcmp(a,'1')
            c=b;
        elseif strcmp(b,'1')
            c=a;
        else
            if ~rise_isa(a,'atom') && ~rise_isa(a,'function')
                a=parenthesize(a,'+-');
            end
            if ~rise_isa(b,'atom') && ~rise_isa(b,'function')
                b=parenthesize(b,'+-');
            end
            [a,b]=commute(a,b);
            c=[a,'*',b];
        end
    case 'power'
        if strcmp(b,'0')||strcmp(a,'1')
            c='1';
        elseif strcmp(b,'1')
            c=a;
        else
            if ~rise_isa(a,'atom') && ~rise_isa(a,'function')
                a=parenthesize(a,'+-*/^');
            end
            if ~rise_isa(b,'atom') && ~rise_isa(b,'function')
                b=parenthesize(b,'+-*/^');
            end
            c=[a,'^',b];
        end
    case {'divide','rdivide'}
        if strcmp(a,'0')
            c='0';
        else
            if ~rise_isa(a,'atom') && ~rise_isa(a,'function')
                a=parenthesize(a,'+-');
            end
            if ~rise_isa(b,'atom') && ~rise_isa(b,'function')
                b=parenthesize(b,'+-*/^');
            end
            c=tryevaluate([a,'/',b]);
        end
    otherwise
        a=tryevaluate(a);
        b=tryevaluate(b);
        c=[func,'(',a,',',b,')'];
end
c=tryevaluate(c);

    function [a,b]=commute(a,b)
        if b(1)<a(1)
            btmp=b;
            atmp=a;
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
end