function rise_startup(flag)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

if nargin<1
    flag=false;
end

% latex and reporting
%--------------------
latex_progs={'pdflatex','epstopdf'};
latex_paths=latex_progs;
retcode=0;
for iprog=1:numel(latex_progs)
    if ispc
        [rcode,latex_paths{iprog}] = system(['findtexmf --file-type=exe ',...
            latex_progs{iprog}]);
    elseif ismac || isunix
        [rcode,latex_paths{iprog}] = ...
            system(['PATH=$PATH:/usr/texbin:/usr/local/bin:/usr/local/sbin;' ...
            'which ',latex_progs{iprog}]);
    else% gnu/linux
        [rcode,latex_paths{iprog}] = system(['which ',latex_progs{iprog}]);
    end
    
    if any(isspace(latex_paths{iprog}))
        latex_paths{iprog}=strcat('"',strtrim(latex_paths{iprog}),'"');
    end
    retcode=retcode||rcode;
end

if retcode
    miktex_version='Not installed (no reporting possible)';
else
    miktex_version=regexp(latex_paths,['(?i)M(iK|ac)TeX(\s+|\',filesep,')(\d+\.\d+)'],'match');
    miktex_version=[miktex_version{:}];
    if isempty(miktex_version)
        miktex_version=regexp(latex_paths,['(?i)M(iK|ac)TeX(\s+|\',filesep,')'],'match');
        miktex_version=[miktex_version{:}];
    end
    if isempty(miktex_version)
        miktex_version='version could not be identified';
    else
        miktex_version=strrep(miktex_version,filesep,' ');
        miktex_version=miktex_version{1};
    end
end
% Decide if using or not the RISE print and plot settings
%--------------------------------------------------------
USE_RISE_PLOT = true;
USE_RISE_PRINT = true;

rise_data=cell(0,2);

if USE_RISE_PRINT
    
    % Version Variables
    % ------------------
    
    NOT_INSTALLED = 'Not installed';
    matlab_version = NOT_INSTALLED;
    optimization_version = [NOT_INSTALLED,'(Fatal error)'];
    statistics_version = NOT_INSTALLED;
    rise_version = NOT_INSTALLED;
    
    vs = ver;
    for jj = 1:length(vs)
        v = vs(jj);
        switch v.Name
            case 'MATLAB'
                matlab_version = [v.Version ' ' v.Release];
            case 'Optimization Toolbox'
                optimization_version = [v.Version ' ' v.Release];
            case {'Statistics Toolbox',...
                    'Statistics and Machine Learning Toolbox'}
                statistics_version = [v.Version ' ' v.Release];
            case 'RISE Toolbox'
                rise_version = [v.Version ' ' v.Release];
        end
    end
    rise_data=[rise_data
        {'matlab_version', matlab_version
        'optimization_version', optimization_version
        'statistics_version', statistics_version
        'rise_version', rise_version
        'rise_required_matlab_version', '7.11'}
        ];
    
    % set page properties for printing
    %---------------------------------
    set(0, 'DefaultFigurePaperOrientation','landscape');
    set(0, 'DefaultFigurePaperType','A4');
    set(0, 'DefaultFigurePaperUnits', 'centimeters');
    set(0, 'DefaultFigurePaperPositionMode', 'manual');
    set(0, 'DefaultFigurePaperPosition', [3.56 2.03 22.56 16.92]);
end


% Plot settings
%--------------
if USE_RISE_PLOT
    
    rise_default_plot_colors={ ...
        [0 0 1],     ...  % 'b'
        [1 0 0],     ...  % 'r'
        [0 1 0],     ...  % 'g'
        [0 0 0],     ...  % 'k'
        [0 1 1],     ...  % 'c'
        [1 0 1],     ...  % 'm'
        [0.565 0.247 0.667],     ...  % pink
        [0.722 0.420 0.274],     ...   % siena
        [0.659 0.541 0.000],     ...   % ocra
        [1 0.604 0.208],     ...   % orange
        [0.502 0.502 0.502],     ...   % dark grey
        [0.733 0.824 0.082],     ...   % ill green
        [0.318 0.557 0.675],     ...   % cobalto
        [0.8 0.2 0.2],     ...
        [0.2 0.2 0.8],     ...
        [0.2 0.9 0.2],     ...
        [0.37 0.9 0.83],   ...
        [0.888 0.163 0.9], ...
        [0 0 0],           ...
        [0 207 255]/255,   ...
        [255 128 0]/255,   ...
        [143 0 0]/255,     ...
        [255 207 0]/255,   ...
        [0.9 0.266 0.593]
        };
    
    rise_data=[rise_data
        {'rise_default_plot_colors',rise_default_plot_colors}];
    set(0, 'DefaultFigureColor', 'w');
end

searchPaths=collect_paths(mfilename);
for jj = 1:numel(searchPaths)
    if flag
        rmpath(searchPaths{jj});
    else
        addpath(searchPaths{jj});
    end
end

%--------------------------------------------------------------------------
target=which(mfilename);
separators=find(target==filesep);
rise_root=target(1:separators(end)-1);
toolbox_folder_name=target(separators(end-1)+1:separators(end)-1);
rise_data=[rise_data
    {'rise_root',rise_root}];
rise_data=[rise_data
    [strcat('rise_',latex_progs(:)),latex_paths(:)]
    {'paths',searchPaths}];

for id=1:size(rise_data,1)
    if flag
        rmappdata(0,rise_data{id,1});
    else
        setappdata(0, rise_data{id,1}, rise_data{id,2});
    end
end

if flag
    set(0,'default');
    % also remove the root path although it was added by the user manually
    %---------------------------------------------------------------------
    rmpath(fileparts(which(mfilename)))
else
    welcome_message();
end

    function welcome_message()
        vv = ver(toolbox_folder_name);
        % paths to documentation should be dynamic
        %-----------------------------------------
        pdf_doc=[rise_root,filesep,'help',filesep,'build',filesep,'latex',filesep,'RISE.pdf'];
        html_doc=[rise_root,filesep,'help',filesep,'build',filesep,'html',filesep,'master_doc.html'];
        
        MiKTeX='MiKTeX';
        if ~ispc
            MiKTeX='MacTeX';
        end
        your_system={
            '* About your system'
            [' - Matlab: ',matlab_version]
            [' - Optimization toolbox: ',optimization_version]
            [' - Statistics toolbox: ',statistics_version]
            [' - ',MiKTeX,': ',miktex_version]
            };
        l1 = '+~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+';
        logo={
            ' _____	  _  ____  _____'
            '|  _  |  (_)|  __||  ___|    |   Rationality In Switching Environments'
            ['| (_) /  | || |__ | |___     |	  Welcome to the ',vv.Name]
            ['|  __ \\  | ||__  ||  ___|    |	  Version: ',vv.Version]
            ['| |  \\ \\ | | __| || |___     |	  Date: ',vv.Date]
            ['|_|   \\_\\|_||____||_____|    |	  Tested with Matlab: ',vv.Release]
            };
        docs={
            ['* Check out the <a href="',strrep(html_doc,'\','\\'),'">html documentation</a>, ',...
            'or the <a href="',strrep(pdf_doc,'\','\\'),'">pdf documentation</a> ']
            };
        concerns={
            ['* Use the <a href="https://github.com/jmaih/RISE_toolbox/issues">',...
            'forum</a> to post concerns, problems, suggestions and desideratas']
            '* Thank you in advance for your feedback !!!'
            };
        tmp=[
            l1
            logo
            your_system
            docs
            concerns
            ];
        if retcode
            tmp=[tmp
                '* pdflatex/epstopdf (Miktex) could not be located!!!'
                ];
        end
        tmp=[tmp
            l1
            ' '];
        try
            fprintf(1,strjoin(tmp(:).','\n'));
        catch
            for irow=1:numel(tmp)
                fprintf(1,[tmp{irow},'\n']);
            end
        end
    end
end

function collect=collect_paths(filename)
fullpath=which(filename);
loc=strfind(fullpath,filename);
fullpath=fullpath(1:loc-2);
subfolders={'m','classes',};
tmp='';
for sfld=1:numel(subfolders)
    tmp=[tmp,genpath([fullpath,filesep,subfolders{sfld}])]; %#ok<AGROW>
end

if ispc
    collect=regexp(tmp,';','split');
elseif ismac || isunix
    collect=regexp(tmp,':','split');
else
    error([mfilename,':: unknown system '])
end

end