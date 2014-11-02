function c=print(deriv,args,long,optimize)
% print - transforms the output of splanar.differentiate into char or functions
%
% Syntax
% -------
% ::
%
%   c=splanar.print(deriv)
%   c=splanar.print(deriv,args)
%   c=splanar.print(deriv,args,long)
%   c=splanar.print(deriv,args,long,optimize)
%
% Inputs
% -------
%
% - **deriv** [output of differentiate]:
%
% - **args** [cellstr|{[]}]: argument list to functions. If empty, print
%   does not create functions. When creating functions, the field
%   'derivatives' is renamed to functions.
%
% - **long** [true|{false}|[]]: if true, functions with shorthands are
%   written with their file names e.g. a+b is written as plus(a,b)
%
% - **optimize** [true|{false}]: optimizes expression by removing e.g.
%   redundant parentheses, etc. This does not work very well so far.
%
% Outputs
% --------
%
% - **c** [structure with same format as deriv]: except that the
%   derivatives field is transformed and the map is populated
%
% More About
% ------------
%
% Examples
% ---------
%
% See also:
if nargin<4
    optimize=false;
    if nargin<3
        long=false;
        if nargin<2
            args=[];
        end
    end
end

nderivs=numel(deriv);

if nderivs>1
    c=deriv(1:0);
    c=update_fieldnames(c,args);
    for id=1:nderivs
        c(id)=splanar.print(deriv(id),args,long,optimize);
    end
    return
end

% char the derivatives
%-----------------------
c=char(deriv.derivatives,long);

% transform to functions if arguments are provided
%--------------------------------------------------
do_functions()

% map the derivatives to their location in the compact matrix
%-------------------------------------------------------------
cmap=[deriv.derivatives.location];
cmapl=cmap(1:2:end);
cmapr=cmap(2:2:end);
cmap=[cmapl(:),cmapr(:)];
% cmap=reshape([deriv.derivatives.location],2,[]);

% format output
%---------------
deriv.derivatives=c;
c=deriv;
c.map=cmap;

% update field names
%--------------------
c=update_fieldnames(c,args);

    function do_functions()
        if ~isempty(args)
            if ischar(args)
                args=cellstr(args);
            end
            args=args(:)';
            
            % put into analytical form
            %--------------------------
            c=parser.analytical_symbolic_form(c,args,'analytic');
            
            % remove unnecessary parentheses
            %-------------------------------
            if optimize
                word='\w+';
                word_par='\w+\(\d+\)';
                c=regexprep(c,['(?<!\w+)(\()(',word,'|',word_par,')(\))'],'$2');
            end
            
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
