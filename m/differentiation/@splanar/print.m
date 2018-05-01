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
%    - **optimize** [true|{false}]: optimizes expression by removing e.g.
%      redundant parentheses, etc. This does not work very well so far.
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
    'deriv',[],@(x)isstruct(x) && all(isfield(x,{'derivatives','nwrt','order','nnz_derivs','partitions'}))
    'args',[],@(x)ischar(x)||iscellstr(x)
    'asfunc', true,@(x)islogical(x)
    'long', false,@(x)islogical(x)
    'optimize', false,@(x)islogical(x)
    };

[deriv,args,asfunc,long,optimize]=...
    parse_arguments(Defaults,'deriv',deriv,'args',args,'asfunc',asfunc,'long',long,'optimize',optimize);

nderivs=numel(deriv);

if nderivs>1
    
    c=deriv(1:0);
    
    c=update_fieldnames(c,args);
    
    for id=1:nderivs
        
        c(id)=splanar.print(deriv(id),args,asfunc,long,optimize);
        
    end
    
    return
    
end

c=[];

cmap=[];

if ~isempty(deriv.derivatives)
    % char the derivatives
    %-----------------------
    c=char(deriv.derivatives,long);
    
    % remove unnecessary parentheses
    %-------------------------------
    do_optimize()
    
    % put into analytic form
    %------------------------
    do_analytic()
    
    % transform to functions if arguments are provided
    %--------------------------------------------------
    c=do_functions(c,args,asfunc);
    
    % map the derivatives to their location in the compact matrix
    %-------------------------------------------------------------
    cmap=[deriv.derivatives.location];
    
    cmapl=cmap(1:2:end);
    
    cmapr=cmap(2:2:end);
    
    cmap=[cmapl(:),cmapr(:)];
    
    % do expansion in case the derivatives are vectorized later on
    %-------------------------------------------------------------
    do_expansion()
    
end
% format output
%---------------
deriv.derivatives=c;

c=deriv;

c.map=cmap;

% update field names
%--------------------
c=update_fieldnames(c,args);

    function do_expansion()
        nd=numel(deriv.derivatives);
        expansions=cell(nd,2); % row, expansion
        offset=0;
        for ider=1:nd
            d=deriv.derivatives(ider);
            row=d.location{1};
            nlocs=numel(d.location{end});
            if d.number_of_columns==1 && nlocs>1
                newguys=offset+ones(nlocs,1);
            else
                newguys=offset+(1:nlocs).';
            end
            expansions(ider,:)={row*ones(nlocs,1),newguys};
            offset=newguys(end);
        end
        expansions=cell2mat(expansions);
        cmapr=cell2mat(cmapr);
        deriv.vectorizer=struct('rows',expansions(:,1),'cols',cmapr(:),'inflator',expansions(:,2));
    end

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
            
            word='\w+';
            
            word_par='\w+\(\d+\)';
            
            c=regexprep(c,['(?<!\w+)(\()(',word,'|',word_par,')(\))'],'$2');
        
        end
        
    end

end

function c=do_functions(c,args,asfunc)

if asfunc && ~isempty(args) && ~isempty(c)
    % this assumes the args are already in correct order/format
    % from the do_analytic part
    
    % build the @() part
    %--------------------
    args=cell2mat(strcat(args,','));
    
    main_string=['@(',args(1:end-1),')'];
    
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

function c=update_fieldnames(c,args)
% reference: http://blogs.mathworks.com/loren/2010/05/13/rename-a-field-in-a-structure-array/
newField='functions';

oldField='derivatives';

if ~isempty(args)
    
    try
        
        [c.(newField)] = c.(oldField);
    
    catch
        % the thing above does not work when c is empty
        try
            
            [c(:).(newField)] = deal(c.(oldField));
        
        catch
            
            [c(:).(newField)] = deal({});
        
        end
        
    end
    
    c = rmfield(c,oldField);

end

end
