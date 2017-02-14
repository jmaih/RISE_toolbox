function disp(t,indent,fullChar)

if nargin < 3
    
    fullChar = false;
    
    if nargin < 2
        
        indent = 4;
        
    end
    
end

between = indent;

within = 2;

% Follow the cmd window's format settings as possible
isLoose = strcmp(get(0,'FormatSpacing'),'loose');

if isLoose
    
    looseline = '\n';
    
else
    
    looseline = '';
    
end

rownames=serial2date(t.date_numbers);

NumberOfVariables=t.NumberOfVariables;

NumberOfObservations=t.NumberOfObservations;

if ~((NumberOfObservations > 0) && (NumberOfVariables > 0))
    
    return
    
end

dblFmt = getFloatFormats();

strongBegin = '';strongEnd = '';

varnameFmt = [strongBegin '%s' strongEnd];

indentSpaces = repmat(' ', NumberOfObservations, indent);   % indent at left margin

betweenSpaces = repmat(' ', NumberOfObservations, between); % space betweeen variables

if isempty(rownames)
    
    tblChars = indentSpaces;
    
else
    
    rownameChars = char(rownames);
    
    rownameWidth = size(rownameChars,2);
    
    rownameChars = strcat(strongBegin,rownameChars,strongEnd);
    
    tblChars = [indentSpaces rownameChars betweenSpaces];
    
end

tblCharsOrig=tblChars;

sizdat=size(t.data);

if numel(sizdat)==2
    
    sizdat=[sizdat,1];
    
end

is_many_pages=sizdat(3)>1;

for ip=1:sizdat(3)
    
    data=mat2cell(t.data(:,:,ip),NumberOfObservations,ones(1,t.NumberOfVariables));
    
    if is_many_pages
        
        fprintf('(:,:,%0.0f) =\n',ip);
        
    end
    
    main_engine(data)
    
end

index(t)

    function main_engine(data)
        
        varnamePads = zeros(1,NumberOfVariables);
        
        for ivar = 1:NumberOfVariables
            
            name = t.varnames{ivar};
            
            var = data{ivar};
            
            if ischar(var)
                
                if ismatrix(var) && (fullChar || (size(var,2) <= 10))
                    % Display individual strings for a char variable that is 2D and no
                    % more than 10 chars.
                    varChars = var;
                    
                else
                    % Otherwise, display a description of the chars.
                    sz = size(var);
                    
                    szStr = ['[1' sprintf('x%d',sz(2:end)) ' char]'];
                    
                    varChars = repmat(szStr,sz(1),1);
                    
                end
                
            else
                % Display the individual data if the var is 2D and no more than 3 columns.
                if ~isempty(var) && ismatrix(var) && (size(var,2) <= 3)
                    
                    if isnumeric(var)
                        
                        varChars = num2str(var,dblFmt);
                        
                    elseif islogical(var)
                        % Display the logical values using meaningful names.
                        strs = ['false'; 'true '];
                        
                        w1 = size(strs,2); w2 = within;
                        
                        varChars = repmat(' ',size(var,1),(size(var,2)-1)*(w1+w2));
                        
                        for j = 1:size(var,2)
                            
                            varChars(:,(j-1)*(w1+w2)+(1:w1)) = strs(var(:,j)+1,:);
                            
                        end
                        
                    else
                        % Display a description of each table element.
                        varChars = getInfoDisplay(var);
                        
                    end
                    
                    % Either the variable is not 2D, or it's empty, or it's too wide
                    % to show. Display a description of each table element.
                else
                    
                    varChars = getInfoDisplay(var);
                    
                end
                
            end
            
            if size(varChars,2) < length(name)
                
                varChars(:,end+1:length(name)) = ' ';
                
            end
            
            varnamePads(ivar) = size(varChars,2)-length(name);
            
                if ivar == 1 % starting over at left margin
                    
                    tblChars = [tblChars varChars]; %#ok<AGROW>
                    
                else
                    
                    tblChars = [tblChars betweenSpaces varChars]; %#ok<AGROW>
                    
                end
                
        end
        
        dispVarNames(true);
        
        disp(tblChars);
        
        dispVarNames();
        
        fprintf(looseline);
        
        if sizdat(3)>1
            
            tblChars=tblCharsOrig;
            
        end
        
        %-----------------------------------------------------------------------
        function dispVarNames(before)
            
            if nargin==0
                
                before=false;
                
            end
            
            if before
                
                do_names()
                
                do_lines()
                
            else
                
                do_lines()
                
                do_names()
                
            end
            
            function do_lines()
                
                if isempty(rownames)
                    
                    fprintf('%s',repmat(' ',1,indent));
                    
                else
                    
                    fprintf('%s',repmat(' ',1,indent+rownameWidth+between));
                    
                end
                
                ii = 1;
                
                ul = repmat('_',1,length(t.varnames{ii})+varnamePads(ii));
                
                fprintf(varnameFmt,ul);
                
                for ii = 2:NumberOfVariables
                    
                    spaces = repmat(' ',1,between);
                    
                    ul = repmat('_',1,length(t.varnames{ii})+varnamePads(ii));
                    
                    fprintf('%s',[spaces sprintf(varnameFmt,ul)]);
                    
                end
                
                fprintf(['\n' looseline]);
                
            end
            
            function do_names()
                
                if isempty(rownames)
                    
                    fprintf('%s',repmat(' ',1,indent));
                    
                else
                    
                    fprintf('%s',repmat(' ',1,indent+rownameWidth+between));
                    
                end
                
                iii = 1;
                
                rightSpaces = repmat(' ',1,ceil(varnamePads(iii) / 2));
                
                leftSpaces = repmat(' ',1,varnamePads(iii)-length(rightSpaces));
                
                fprintf(varnameFmt,[leftSpaces t.varnames{iii} rightSpaces]);
                
                for iii = 2:NumberOfVariables
                    
                    rightSpaces = repmat(' ',1,ceil(varnamePads(iii) / 2));
                    
                    leftSpaces = repmat(' ',1,varnamePads(iii)-length(rightSpaces)+between);
                    
                    fprintf(varnameFmt,[leftSpaces t.varnames{iii} rightSpaces]);
                    
                end
                
                fprintf('\n');
                
            end
            
        end
        
    end % main engine

end % main function

%-----------------------------------------------------------------------
function [dblFmt] = getFloatFormats()
% Display for double/single will follow 'format long/short g/e' or 'format bank'
% from the command window. 'format long/short' (no 'g/e') is not supported
% because it often needs to print a leading scale factor.
switch get(0,'Format')
    
    case {'short' 'shortG' 'shortEng'}
        
        dblFmt  = '%.5g    ';
        
    case {'long' 'longG' 'longEng'}
        
        dblFmt  = '%.15g    ';
        
    case 'shortE'
        
        dblFmt  = '%.4e    ';
        
    case 'longE'
        
        dblFmt  = '%.14e    ';
        
    case 'bank'
        
        dblFmt  = '%.2f    ';
        
    otherwise % rat, hex, + fall back to shortg
        
        dblFmt  = '%.5g    ';
        
end

end

%-----------------------------------------------------------------------
% function str = removeBraces(str)
% str = regexprep(str,'\{(.*)\}','$1');
% end

%-----------------------------------------------------------------------
function varChars = getInfoDisplay(var)

sz = size(var);

if ismatrix(var)
    
    szStr = ['[1' sprintf('x%d',sz(2:end))];
    
else
    
    szStr = ['[1' sprintf('x%d',sz(2)) sprintf('x%d',sz(3:end))];
    
end

varChars = repmat([szStr ' ' class(var) ']'],sz(1),1);

end
