function table_displayer(the_data,colnames,rownames,prologue,epilogue,indent,fullChar)

if nargin < 7
    
    fullChar = false;
    
    if nargin < 6
        
        indent = 4;
        
        if nargin <5
            
            epilogue=[];
            
            if nargin<4
                
                prologue=[];
                
            end
            
        end
        
    end
    
end

cell_style=iscell(the_data);

if cell_style
    
    NumberOfObservations=size(the_data{1},1);
    
    NumberOfVariables=numel(the_data);
    
    NumberOfPages=1;
    
else
    
[NumberOfObservations,NumberOfVariables,NumberOfPages]=size(the_data);

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

is_many_pages=NumberOfPages>1;

for iiii=1:numel(prologue)
    
    fprintf(1,'%s \n',prologue{iiii});
    
end

for ip=1:NumberOfPages
    
    if cell_style
        
        data=the_data;
        
    else
        
        data=mat2cell(the_data(:,:,ip),NumberOfObservations,...
            ones(1,NumberOfVariables));
        
    end
    
    if is_many_pages
        
        fprintf('(:,:,%0.0f) =\n',ip);
        
    end
    
    main_engine(data)
    
end

for iiii=1:numel(epilogue)
    
    fprintf(1,'%s \n',epilogue{iiii});
    
end

    function main_engine(data)
        
        varnamePads = zeros(1,NumberOfVariables);
        
        for ivar = 1:NumberOfVariables
            
            name = colnames{ivar};
            
            datai = data{ivar};
            
            if ischar(datai)
                
                if ismatrix(datai) && (fullChar || (size(datai,2) <= 10))
                    % Display individual strings for a char variable that is 2D and no
                    % more than 10 chars.
                    dataChars = datai;
                    
                else
                    % Otherwise, display a description of the chars.
                    sz = size(datai);
                    
                    szStr = ['[1' sprintf('x%d',sz(2:end)) ' char]'];
                    
                    dataChars = repmat(szStr,sz(1),1);
                    
                end
                
            else
                % Display the individual data if the var is 2D and no more than 3 columns.
                if ~isempty(datai) && ismatrix(datai) && (size(datai,2) <= 3)
                    
                    if isnumeric(datai)
                        
                        dataChars = num2str(datai,dblFmt);
                        
                    elseif islogical(datai)
                        % Display the logical values using meaningful names.
                        strs = ['false'; 'true '];
                        
                        w1 = size(strs,2); w2 = within;
                        
                        dataChars = repmat(' ',size(datai,1),(size(datai,2)-1)*(w1+w2));
                        
                        for j = 1:size(datai,2)
                            
                            dataChars(:,(j-1)*(w1+w2)+(1:w1)) = strs(datai(:,j)+1,:);
                            
                        end
                        
                    else
                        % Display a description of each table element.
                        dataChars = getInfoDisplay(datai);
                        
                    end
                    
                    % Either the variable is not 2D, or it's empty, or it's too wide
                    % to show. Display a description of each table element.
                else
                    
                    dataChars = getInfoDisplay(datai);
                    
                end
                
            end
            
            if size(dataChars,2) < length(name)
                
                dataChars(:,end+1:length(name)) = ' ';
                
            end
            
            varnamePads(ivar) = size(dataChars,2)-length(name);
            
                if ivar == 1 % starting over at left margin
                    
                    tblChars = [tblChars dataChars]; %#ok<AGROW>
                    
                else
                    
                    tblChars = [tblChars betweenSpaces dataChars]; %#ok<AGROW>
                    
                end
                
        end
        
        dispVarNames(true);
        
        disp(tblChars);
        
        dispVarNames();
        
        fprintf(looseline);
        
        if NumberOfPages>1
            
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
                
                ul = repmat('_',1,length(colnames{ii})+varnamePads(ii));
                
                fprintf(varnameFmt,ul);
                
                for ii = 2:NumberOfVariables
                    
                    spaces = repmat(' ',1,between);
                    
                    ul = repmat('_',1,length(colnames{ii})+varnamePads(ii));
                    
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
                
                fprintf(varnameFmt,[leftSpaces colnames{iii} rightSpaces]);
                
                for iii = 2:NumberOfVariables
                    
                    rightSpaces = repmat(' ',1,ceil(varnamePads(iii) / 2));
                    
                    leftSpaces = repmat(' ',1,varnamePads(iii)-length(rightSpaces)+between);
                    
                    fprintf(varnameFmt,[leftSpaces colnames{iii} rightSpaces]);
                    
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