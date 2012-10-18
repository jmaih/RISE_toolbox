function string=char(obj,unravel,isparent)
if nargin<3
    isparent=false;
    if nargin<2
        unravel=false;
    end
end
% char itself is already taken care of
args=reprocess_arguments(obj.args);
if isa(obj.name,'double')
    string=num2str(obj.name,10);
elseif isempty(args) % variable
    string=obj.name; % this should be a char
elseif ~unravel && ~isempty(obj.ref) && ~isparent
    string=obj.ref;
else
    switch obj.name
        case {'plus','minus','times','power'}
            operator=str2func(['my',obj.name]);
            string=operator(mychar(args{1},unravel),mychar(args{2},unravel));
        case 'uplus'
            string=myplus('0',mychar(args{1},unravel));
        case 'uminus'
            string=myminus('0',mychar(args{1},unravel));
        case {'mtimes'}
            string=mytimes(mychar(args{1},unravel),mychar(args{2},unravel));
        case {'mpower'}
            string=mypower(mychar(args{1},unravel),mychar(args{2},unravel));
        case {'rdivide','mrdivide'}
            string=mydivide(mychar(args{1},unravel),mychar(args{2},unravel));
        case {'min','max','gt','lt','ge','le'}
            string=[obj.name,'(',mychar(args{1},unravel),',',mychar(args{2},unravel),')'];
        case {'ldivide','mldivide'}
            string=mydivide(mychar(args{2},unravel),mychar(args{1},unravel));
        case {'exp','log','log10','sin','asin','sinh','cos','acos','cosh',...
                'tan','atan','tanh','abs','sqrt','isreal','sign'}
            string=[obj.name,'(',mychar(args{1},unravel),')'];
        case {'normpdf','normcdf'}
            string=[obj.name,'(',mychar(args{1},unravel),',',mychar(args{2},unravel),',',mychar(args{3},unravel),')'];
    end
end
if isparent
    string=[obj.ref,'=',string];
end
end

function args=reprocess_arguments(args)
for ii=1:numel(args)
    if isnumeric(args{ii})
        args{ii}=num2str(args{ii},10);
    end
end
end

function cc=mychar(obj,unravel)
if ischar(obj)
    cc=obj;
else
    cc=char(obj,unravel);
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
flag=isempty(regexp(string,'[/*\-+^]','start'));
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

