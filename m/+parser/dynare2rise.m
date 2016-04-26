function dynare2rise(dynFileName,riseFileName,stderr_name,fast)

if nargin<4
    
    fast=[];
    
    if nargin<3
        
        stderr_name=[];
        
    end
    
end

if isempty(fast)
    
    fast=false;
    
end

if isempty(stderr_name)
    
    stderr_name='std';
    
end

dynFileName = strtrim(dynFileName);

riseFileName = strtrim(riseFileName);

% read input file
%-----------------
raw_code = read_file();

% remove block comments
%-----------------------
rise_code = regexprep(raw_code,'/\*(.*)\*/','');

header = sprintf('Conversion of Dynare file [%s] into RISE file [%s]\n',...
    dynFileName,riseFileName);

% Remove line comments.
rise_code = regexprep(rise_code,'//(.*?\n)','');

% Convert char(10) to white space.
eol=char(10);
rise_code = converteols(rise_code);

% remove multiple white characters
rise_code = regexprep(rise_code,'\s+',' ');

rise_code = regexprep(rise_code,'model\(linear\)','model');

% predetermined variables
%-------------------------
pred_vars=get_list('predetermined_variables',';');

% endogenous: var section
%-------------------------
endo_list=get_list('var',';');

% exogenous: varexo section
%----------------------------
exo_list=get_list('varexo',';');

% parameters: parameters section
%--------------------------------
par_list=get_list('parameters',';');

% add new parameters
par_list=[par_list,strcat([stderr_name,'_'],exo_list)];

% model equations
%-----------------
model_eqtns = get_list('model;','end;',false);

% replace time indices
model_eqtns = regexprep(model_eqtns,'(?<=\w)\(([\+-]?\d*)\)','{$1}');

% take care of predetermined variables
do_predetermined();

% steady state model equations
%------------------------------
ssmodel_eqtns = get_list('steady_state_model;','end;',false);

% push standard deviations into equations directly
%-------------------------------------------------
add_stdev_to_model(fast)

% Looking for covariances: shocks section
%----------------------------------------
shocks_block = get_list('shocks;','end;',false);

if ~isempty(shocks_block)
    
    match = regexpi(shocks_block,'var\s*\w+\s*,\s*\w+\s*=.*?;','match');
    
    for ii = 1 : length(match)
        
        header = [header,sprintf('Cross-covariance ignored:\n  %s\n',...
            match{ii})]; %#ok<AGROW>
        
    end
    
end

% Create and save RISE code.
%---------------------------
timestamp = datestr(now);

recreate_code();

    function do_predetermined()
        
        if ~isempty(pred_vars)
            
            just_neat=@myreplace; %#ok<NASGU>
            
            % add {0} to guys that are current
            %----------------------------------
            xpress0=parser.cell2matize(pred_vars);
            xpress=['(?<!\w+)',xpress0,'(?!\w+)(?!{)'];
            model_eqtns=regexprep(model_eqtns,xpress,'$1{0}');
            
            % substract 1 to the time indices
            %---------------------------------
            xpress=['(?<!\w+)',xpress0,'(?!\w+){(\+|\-)?(\d+)}'];
            repl='${just_neat($1,$2,$3)}';
            model_eqtns=regexprep(model_eqtns,xpress,repl);
            
            % remove {0}
            %------------
            model_eqtns=strrep(model_eqtns,'{0}','');
                        
        end
        
        function b=myreplace(a1,a2,a3)
            
            b=a1;
            
            a2a3_1=str2double([a2,a3])-1;
            
            if a2a3_1
                
                b=[b,'{',int2str(a2a3_1),'}'];
                
            end
            
        end
        
    end

    function add_stdev_to_model(fast)
        
        if fast
            
            for ishock=1:numel(exo_list)
                
                expre=exo_list{ishock};
                
                repl=[stderr_name,'_',expre,'*',expre];
                
                model_eqtns=regexprep(model_eqtns,expre,repl);
                
            end
            
        else
            
            expre=cell2mat(strcat(exo_list,'|'));
            
            expre=['(?<!\w+)(',expre(1:end-1),')(?!\w+)'];
            
            repl=[stderr_name,'_$1*$1'];
            
            model_eqtns=regexprep(model_eqtns,expre,repl);
        end
    end

    function rise_code = recreate_code()
        
        header = [header,sprintf('\n Done %s.',timestamp)];
        
        rise_code = [
            'endogenous ',endo_list,eol,...
            'exogenous ',exo_list,eol,...
            'parameters ',par_list,eol,...
            'model ',model_eqtns,eol
            ];
        
        if ~isempty(ssmodel_eqtns)
            rise_code = [rise_code,...
                'steady_state_model ',ssmodel_eqtns,eol
                ];
        end
        
        rise_code = regexprep(rise_code,'(\n)\s*','$1  ');
        
        rise_code = regexprep(rise_code,';',sprintf(';\n\n'));
        
        rise_code = [
            regexprep(['%',header],'(\n)','$1%'),...
            rise_code
            ];
        
        write2file(rise_code,riseFileName);
        
    end

    function list=get_list(typeof,closing,do_match)
        
        if nargin<3
            
            do_match=true;
            
        end
        
        tokens = regexpi(rise_code,...
            [typeof,'\s*(.*?)',closing],'tokens','once');
        
        if isempty(tokens)
            
            list=tokens;
            
            return
            
        end
        
        if do_match
            
            list = regexp(tokens{1},'\w*','match');
            
        else
            
            list = tokens{1};
            
        end
        
    end

    function raw_code = read_file()
        
        fid = fopen(dynFileName,'r');
        
        if fid == -1
            
            if exist(dynFileName,'file') == false
                
                error('Unable to find ''%s''.',dynFileName);
                
            else
                
                error('Unable to open ''%s'' for reading.',dynFileName);
                
            end
            
        end
        
        raw_code = transpose(fread(fid,'*char'));
        
        fclose(fid);
        
    end

end

function x = converteols(x)
%x = regexprep(x,'\r\n?','\n');
% This is much faster:
% Windows:
x = strrep(x,sprintf('\r\n'),sprintf('\n'));
% Apple:
x = strrep(x,sprintf('\r'),sprintf('\n'));

end

function write2file(C,filename)

fid = fopen(filename,'w+');

if fid == -1
    
    error(['Cannot open file ''%s'' for writing.',filename]);
    
end

if iscellstr(C)
    
    C = sprintf('%s\n',C{:});
    
    if ~isempty(C)
        
        C(end) = '';
        
    end
    
end

count = fwrite(fid,C,'char');

if count ~= length(C)
    
    fclose(fid);
    
    error(['Cannot write character string to file ''%s''.',filename]);
    
end

fclose(fid);

end