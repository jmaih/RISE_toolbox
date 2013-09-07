function c=char(x,long)
if nargin<2
    long=false;
end

if isnumeric(x)
    n=numel(x.func);
    c=num2str(x.func);
    if n>1
        c=['[',c,']'];
    end
elseif isempty(x.args)
    c=x.func;
else
    if ~any(x.incidence)
        c='0';
        return
    end
    args=x.args;
    for iarg=1:numel(args)
        args{iarg}=char(args{iarg},long);
    end
    switch x.func
        % unary functions
        %----------------
        case {'abs','acos','acosh','asin','asinh','atan','atanh','cos',...
                'cosh','cot','erf','exp','log','log10','sign','sin','sinh',...
                'tan'}
            c=[x.func,'(',args{1},')'];
        case 'uplus'
            c=args{1};
        case 'uminus'
            if long,c=['uminus(',args{1},')']; else c=['-(',args{1},')'];end
            % binary functions
            %-----------------
        case 'and'
            if long,c=['and(',args{1},',',args{2},')'];else c=['(',args{1},' & ',args{2},')'];end
        case 'eq'
              if long,c=['eq(',args{1},',',args{2},')'];else c=['(',args{1},'==',args{2},')'];end
        case 'ge'
              if long,c=['ge(',args{1},',',args{2},')'];else c=['(',args{1},'>=',args{2},')'];end
        case 'gt'
               if long,c=['gt(',args{1},',',args{2},')'];else c=['(',args{1},'>',args{2},')'];end
       case {'if_elseif','if_then_else'}
            c=[x.func,'('];
            for iarg=1:numel(args)
                c=[c,args{iarg},','];
            end
            c=[c(1:end-1),')'];
        case 'le'
              if long,c=[x.func,'(',args{1},',',args{2},')'];else c=['(',args{1},'<=',args{2},')'];end
        case 'lt'
              if long,c=[x.func,'(',args{1},',',args{2},')'];else c=['(',args{1},'<',args{2},')'];end
        case {'max','min'}
            c=[x.func,'(',args{1},',',args{2},')'];
        case 'minus'
              if long,c=[x.func,'(',args{1},',',args{2},')'];else c=[args{1},'-(',args{2},')'];end
        case 'mpower'
            % use the vectorized form directly
              if long,c=['power(',args{1},',',args{2},')'];else c=['(',args{1},').^(',args{2},')'];end
        case 'mrdivide'
            % use the vectorized form directly
              if long,c=['divide(',args{1},',',args{2},')'];else c=['(',args{1},')./(',args{2},')'];end
        case 'mtimes'
            % use the vectorized form directly
              if long,c=['times(',args{1},',',args{2},')'];else c=['(',args{1},').*(',args{2},')'];end
        case 'ne'
               if long,c=[x.func,'(',args{1},',',args{2},')'];else c=['(',args{1},'~=',args{2},')'];end
       case 'or'
              if long,c=[x.func,'(',args{1},',',args{2},')'];else c=['(',args{1},'|',args{2},')'];end
        case 'plus'
              c=[args{1},'+',args{2}];
        case {'normalcdf','normalpdf'}
              c=[x.func,'(',args{1},',',args{2},',',args{3},')'];
        otherwise
    end
end