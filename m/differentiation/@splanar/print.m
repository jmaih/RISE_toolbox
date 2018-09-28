function c=print(deriv,args,asfunc,long,optimize)
% print - transforms the output of splanar.differentiate into char or functions
%
% ::
%
%
%   c=splanar.print(deriv)
%   c=splanar.print(deriv,args)
%   c=splanar.print(deriv,args,asfunc)
%   c=splanar.print(deriv,args,asfunc,long)
%   c=splanar.print(deriv,args,asfunc,long,optimize)
%
% Args:
%
%    - **deriv** [output of differentiate]:
%
%    - **args** [cellstr|{[]}]: argument list to functions. If empty, print
%      does not create functions. When creating functions, the field
%      'derivatives' is renamed to functions.
%
%    - **long** [true|{false}|[]]: if true, functions with shorthands are
%      written with their file names e.g. a+b is written as plus(a,b)
%
%    - **asfunc** [{true}|false|[]]: if true, the derivatives are put into
%      anonymous functions. Else they are just returned as char or cellstr.
%
%    - **optimize** [{true}|false]: optimizes expression by removing e.g.
%      redundant parentheses, etc.
%
% Returns:
%    :
%
%    - **c** [structure with same format as deriv]: except that the
%      derivatives field is transformed and the map is populated
%
% Note:
%
% Example:
%
%    See also:

if nargin<5
    
    optimize=[];
    
    if nargin<4
        
        long=[];
        
        if nargin<3
            
            asfunc=[];
            
            if nargin<2
                
                args=[];
                
            end
            
        end
        
    end
    
end

Defaults={
    'deriv',[],@(x)isstruct(x) && all(isfield(x,{'derivatives','nwrt',...
    'order','nnz_derivs'}))
    'args',[],@(x)ischar(x)||iscellstr(x)
    'asfunc', true,@(x)islogical(x)
    'long', false,@(x)islogical(x)
    'optimize', true,@(x)islogical(x)
    };

[deriv,args,asfunc,long,optimize]=...
    parse_arguments(Defaults,'deriv',deriv,'args',args,'asfunc',asfunc,...
    'long',long,'optimize',optimize);

nderivs=numel(deriv);

if nderivs>1
    
    c=deriv(1:0);
    
    for id=1:nderivs
        
        c(id)=splanar.print(deriv(id),args,asfunc,long,optimize);
        
    end
    
    return
    
end

c=[];

% build the @() part
%--------------------
main_string=[];

if ~isempty(args)
    
    args_=cell2mat(strcat(args,','));
    
    main_string=['@(',args_(1:end-1),')'];
    
end

if ~isempty(deriv.derivatives)
    % char the derivatives
    %-----------------------
    c=char(deriv.derivatives,long,main_string);
    
    % remove unnecessary parentheses
    %-------------------------------
    do_optimize()
    
    % put into analytic form
    %------------------------
    do_analytic()
    
    % transform to functions if arguments are provided
    %--------------------------------------------------
    c=do_functions(c,main_string,asfunc);
    
    % map the derivatives to their location in the extended matrix
    %-------------------------------------------------------------
    if ~iscell(c)
        
        c={c};
        
    end
    
    c=[c,c];
    
    for ii=1:size(c,1)
        
        c{ii,2}=deriv.derivatives(ii).location;
        
    end
    
end
% format output
%---------------
deriv.derivatives=c;

c=deriv;

    function do_analytic()
        
        if ~isempty(args) && ~isempty(c)
            
            if ischar(args)
                
                args=cellstr(args);
                
            end
            
            args=args(:)';
            
            % put into analytical form
            %--------------------------
            c=parser.analytical_symbolic_form(c,args,'analytic');
            
        end
        
    end

    function do_optimize()
        
        if optimize
            
            pattern=['\<',... beginning of word
                '(\()',... capture (
                '(',... beginning of group
                '\[[0-9,\-.]+\]|',... vectors of numbers
                '\w+(\([\w.\-\+\^/]+\))?',... words possibly followed by (
                ')',... end of group
                '(\))'];
            
            c=regexprep(c,pattern,'$2');
            
        end
        
    end

end

function c=do_functions(c,main_string,asfunc)

if asfunc && ~isempty(main_string) && ~isempty(c)
    % this assumes the args are already in correct order/format
    % from the do_analytic part
    
    % concatenante
    %--------------
    c=strcat(main_string,c);
    
    char_type=ischar(c);
    
    if char_type
        
        c={c};
        
    end
    
    % build the functions
    %--------------------
    c=cellfun(@(x)str2func(x),c,'uniformOutput',false);
    
    if char_type
        
        c=c{1};
        
    end
    
end

end
