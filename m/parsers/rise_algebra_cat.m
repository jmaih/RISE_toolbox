function c=rise_algebra_cat(fun,a,b,debug)
if nargin<4
    debug=false;
    if nargin<3
        b=[];
        if nargin<2
            error('number of arguments should be at least 2')
        end
    end
end
if ischar(a)
    a=cellstr(a);
end
if ~isempty(b)||ischar(b)
    b=cellstr(b);
end
switch fun
    case 'uplus'
        c=a;
        return
    case 'uminus'
        sab=numel(a);
        fun='-';
        [aa,~,~,~,same_a]=format_argument(a,sab);
        c=strcat('-',parenth(aa,same_a));
        return
    case 'plus'
        fun='+';
    case 'minus'
        fun='-';
    case {'mtimes','times'}
        fun='*';
    case {'mrdivide','rdivide'}
        fun='/';
    case {'mpower','power'}
        fun='^';
    case {'+','-','*','/','^'};
        % do nothing
    otherwise
        error('unknown type of binary function')
end

sa=size(a);
sb=size(b);
if sa(1)>1 || sb(1)>1
    error('arguments to concatenate should be ROW vectors. Use loop if necessary')
end
sab=max(sa(2),sb(2));
if ~any(sa(1)==[1 sab])||~any(sb(1)==[1 sab])
    error('vectors can only have one column or the same number of columns')
end
c={'0'};
c=c(1,ones(1,sab));
[aa,ar_l,ar_0,ar_1,same_a]=format_argument(a,sa(2));
[bb,br_l,br_0,br_1,same_b]=format_argument(b,sb(2));
ab_r=strcmp(aa,bb);

no_letter=~(ar_l|br_l);
c(no_letter)=myeval(strcat('(',aa(no_letter),')',fun,'(',bb(no_letter),')'));

if all(no_letter)
    return
end

% remove no_letter location as it has already been processed
ar_0(no_letter)=false; ar_1(no_letter)=false;
br_0(no_letter)=false; br_1(no_letter)=false;
ab_r(no_letter)=false;
switch fun
    case '+'
        myplus()
    case '-'
        myminus()
    case '*'
        mytimes()
    case '/'
        mydivide()
    case '^'
        mypower()
    otherwise
        error(['function ',fun,' is undefined'])
end

    function [x,x_l,x_0,x_1,same]=format_argument(x,ncols)
        x_l=myisletter(x);
        x_0=strcmp(x,'0');
        x_1=strcmp(x,'1');
        if ncols<sab
            x=x(1,ones(1,sab));
            x_0=x_0(1,ones(1,sab));
            x_1=x_1(1,ones(1,sab));
            x_l=x_l(1,ones(1,sab));
            same=true;
        else
            same=all(strcmp(x{1},x));
        end
    end

    function d=parenth(d,same_shit) 
        dejavu=cell(2,0);
        [dejavu,d(1)]=myparenth(dejavu,d(1));
        for kk=2:numel(d)
            if same_shit
                if kk==2
                    d(kk:end)=d(1);
                end
            else
                [dejavu,d(kk)]=myparenth(dejavu,d(kk));
            end
        end
        function [dejavu,dd]=myparenth(dejavu,dd)
            loc=find(strcmp(dd,dejavu(1,:)));
            if isempty(loc)
                old_d=dd;
                dd=parenthesize(dd,fun,debug);
                dejavu=[dejavu,{old_d;dd}];
            else
                dd=dejavu(2,loc);
            end
        end
    end

    function mypower()
        c(br_0)={'1'};
        br_1(br_0)=false;
        c(br_1)=aa(br_1);
        rest=~(no_letter|br_0|br_1);
        c(rest)=strcat(parenth(aa(rest),same_a),'^',parenth(bb(rest),same_b));
    end

    function mydivide()
        br_1(ar_0)=false;
        ab_r(ar_0)=false;
        c(br_1)=aa(br_1);
        ab_r(br_1)=false;
        c(ab_r)={'1'};
        rest=~(no_letter|ar_0|br_1|ab_r);
        c(rest)=strcat(parenth(aa(rest),same_a),'/',parenth(bb(rest),same_b));
    end

    function mytimes()
        c(ar_1)=bb(ar_1);
        c(br_1)=aa(br_1);
        ab_r(br_1)=false;
        c(ab_r)=strcat('(',aa(ab_r),')','^2');
        rest=~(no_letter|ab_r|ar_1|br_1|ar_0|br_0);
        if any(rest)
            a_rest=parenth(aa(rest),same_a);
            b_rest=parenth(bb(rest),same_b);
            c(rest)=strcat(a_rest,'*',b_rest);
        end
    end

    function myminus()
        c(ab_r)={'0'};
        br_0(ab_r)=false;
        c(br_0)=aa(br_0);
        rest=~(no_letter|ab_r|br_0);
        if any(rest)
            b_rest=parenth(bb(rest),same_b);
            a_rest=parenth(aa(rest),same_b);
            arz=strcmp('0',a_rest);
            c(rest(arz))=strcat('-',b_rest(arz));
            c(rest(~arz))=strcat(a_rest(~arz),'-',b_rest(~arz));
        end
    end

    function myplus()
        c(ar_0)=bb(ar_0);
        br_0(ar_0)=false;
        ab_r(ar_0)=false;
        c(br_0)=aa(br_0);
        ab_r(br_0)=false;
        c(ab_r)=strcat('2*(',aa(ab_r),')');
        rest=~(no_letter|ar_0|br_0|ab_r);
        c(rest)=strcat(aa(rest),'+',bb(rest));
    end
end

function flag=myisletter(y)
flag=cellfun(@(x)any(isletter(x)),y,'UniformOutput',true);
end

function c=myeval(y)
% push the results into cell arrays directly
c=cellfun(@(x)sprintf('%0.10g',eval(x)),y,'UniformOutput',false);
end
