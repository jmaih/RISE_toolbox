function c=char(x,long)
% char - transforms splanar derivatives to strings
%
% ::
%
%
%   c=char(x)
%   c=char(x,long)
%
% Args:
%
%    - **x** [splanar|cell array of splanar] : scalar or vector of splanar
%      objects or cell array of splanar objects
%
%    - **long** [true|{false}] : if true, all functions are written with their
%      function names and not with their shorthand. e.g. a+b will be written
%      plus(a,b)
%
% Returns:
%    :
%
%    - **c** [char|cellstr]: derivatives
%
% Note:
%
% Example:
%
%    See also:

if nargin<2
    
    long=false;

end

if iscell(x)
    
    x=[x{:}];

end

if ~isa(x,'splanar')
    
    error('first input argument must be a splanar or a cell array of splanar objects')

end

nx=numel(x);

if nx>1
    
    c=cell(nx,1);
    
    for irow=1:nx
        
        c{irow}=char(x(irow),long);
    
    end
    
    return

end

if isnumeric(x)
    
    n=numel(x.func);
    
    % c=num2str(x.func);
    c=sprintf('%0.15g',full(x.func(1)));
    
    for ii=2:n
        
        c=[c,',',sprintf('%0.15g',full(x.func(ii)))]; %#ok<AGROW>
    
    end
    
    if n>1
        
        c=['[',c,']'];
    
    end
    
elseif isempty(x.args)
    
    c=x.func;

else
    
    args=x.args;
   
    for iarg=1:numel(args)
        
        args{iarg}=char(args{iarg},long);
    
    end
    
    the_func=x.func;
    % use vectorized operation directly
    %----------------------------------
    if any(strcmp(the_func,{'mpower','mrdivide','mtimes'}))% <-- ismember(the_func,{'mpower','mrdivide','mtimes'})
        
        the_func=the_func(2:end);
    
    end
    
    if strcmp(the_func,'plus')
        
        c=[args{1},'+',args{2}];
    
    elseif long||any(strcmp(the_func,{'abs','acos','acosh','asin','asinh','atan','atanh','cos',...
            'cosh','cot','erf','exp','log','log10','sign','sin','sinh',...
            'tan','if_elseif','if_then_else','normcdf','normpdf','norminv',...
            'max','min','sqrt','steady_state','betapdf','betacdf',...
            'betainv'}))
        
        c=cell2mat(strcat(args,','));
        
        c=[the_func,'(',c(1:end-1),')'];
    
    else
        
        switch the_func
            
            case 'uplus'
                
                c=args{1};
            
            case 'uminus'
                
                c=['-(',args{1},')'];
            
            case 'and'
                
                c=['(',args{1},' & ',args{2},')'];
            
            case 'eq'
                
                c=['(',args{1},')==(',args{2},')'];
            
            case 'ge'
                
                c=['(',args{1},')>=(',args{2},')'];
            
            case 'gt'
                
                c=['(',args{1},')>(',args{2},')'];
            
            case 'le'
                
                c=['(',args{1},')<=(',args{2},')'];
            
            case 'lt'
                
                c=['(',args{1},')<(',args{2},')'];
            
            case 'minus'
                
                c=[args{1},'-(',args{2},')'];
            
            case 'power'
                % use the vectorized form directly
                c=['(',args{1},').^(',args{2},')'];
            
            case 'rdivide'
                % use the vectorized form directly
                c=['(',args{1},')./(',args{2},')'];
            
            case 'times'
                % use the vectorized form directly
                c=['(',args{1},').*(',args{2},')'];
            
            case 'ne'
                
                c=['(',args{1},')~=(',args{2},')'];
            
            case 'or'
                
                c=['(',args{1},')|(',args{2},')'];
            
            otherwise
                
                error(['unknown function :: ',the_func])
        
        end
        
    end
    
end

end