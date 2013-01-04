function str=neat(func,varargin)

nargs=length(varargin);
if nargs==0
    str=func;
else
    varargin{1}=varargin{1};
    if nargs>1
        v=varargin{2};
        if nargs>2
            w=varargin{3};
        end
    end
    really_neat=true;
    
    switch func
        case {'+','plus'}
            if (strcmp(varargin{1},'0') && strcmp(v,'0'))
                str='0';
            elseif strcmp(varargin{1},'0')
                str=v;
            elseif strcmp(varargin{1},v)
                if really_neat && check_offend(varargin{1},'+-')
                    varargin{1}=['(',varargin{1},')'];
                end
                str=['2*',varargin{1}];
            elseif strcmp(v,'0')
                str=varargin{1};
            else
                str=[varargin{1},'+',v];
            end
        case {'-','minus'}
            if strcmp(varargin{1},v)
                str='0';
            elseif strcmp(varargin{1},'0')
                if really_neat && check_offend(v,'+-')
                    v=['(',v,')'];
                end
                str=['-',v];
            elseif strcmp(v,'0')
                str=varargin{1};
            else
                if really_neat && check_offend(v,'+-')
                    v=['(',v,')'];
                end
                str=[varargin{1},'-',v];
            end
        case {'*','mtimes','times'}
            if strcmp(varargin{1},'0') || strcmp(v,'0')
                str='0';
            elseif strcmp(varargin{1},'1') && strcmp(v,'1')
                str='1';
            elseif strcmp(varargin{1},'1')
                str=v;
            elseif strcmp(v,'1')
                str=varargin{1};
            elseif strcmp(varargin{1},v)
                if really_neat  && check_offend(varargin{1},'+-/*^')
                    varargin{1}=['(',varargin{1},')'];
                end
                str=[varargin{1},'.^2'];
            else
                if really_neat
                    if check_offend(varargin{1},'+-')
                        varargin{1}=['(',varargin{1},')'];
                    end
                    if check_offend(v,'+-')
                        v=['(',v,')'];
                    end
                end
                str=[varargin{1},'.*',v];
            end
        case {'/','mrdivide','rdivide'}
            if strcmp(varargin{1},'0')
                str='0';
            elseif strcmp(varargin{1},v)
                str='1';
            elseif strcmp(v,'1')
                str=varargin{1};
            else
                if really_neat
                    if check_offend(varargin{1},'+-')
                        varargin{1}=['(',varargin{1},')'];
                    end
                    if check_offend(v,'+-*/')
                        v=['(',v,')'];
                    end
                end
                str=[varargin{1},'./',v];
            end
        case {'^','mpower','power'}
            if strcmp(v,'0')
                str='1';
            elseif strcmp(v,'1')||strcmp(varargin{1},'1')
                str=varargin{1};
            else
                if really_neat
                    if check_offend(varargin{1},'+-*/^')
                        varargin{1}=['(',varargin{1},')'];
                    end
                    if check_offend(v,'+-*/^')
                        v=['(',v,')'];
                    end
                end
                str=[varargin{1},'.^',v];
            end
        case 'uplus'
            str=varargin{1};
        case 'uminus'
            if really_neat
                if check_offend(varargin{1},'+-')
                    varargin{1}=['(',varargin{1},')'];
                end
            end
            str=['-',varargin{1}];
        case {'exp','log','cos','acos','cosh','sin','asin','sinh','tan',...
                'atan','tanh','sqrt'}
            str=[func,'(',varargin{1},')'];
        case {'normpdf','normcdf'}
            str=[func,'(',varargin{1},',',v,',',w,')'];
        otherwise
%             str=[func,'(',varargin{:},')'];
            error([func,' is undefined for objects of class ',mfilename])
    end
end
end

function flag=check_offend(v,offend)
flag=false;
iter=0;
while iter<length(v)
    iter=iter+1;
    if any(v(iter)==offend)
        flag=true;
        break
    end
end
% flag=any(ismember(v,offend));
end


