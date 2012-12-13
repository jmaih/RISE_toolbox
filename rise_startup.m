function rise_startup(flag)
if nargin<1
    flag=false;
end

% this function sets up RISE (it replaces setpaths)

[retcode,pdflatex_path] = system('findtexmf --file-type=exe pdflatex');

%-----------------------------------------------------------------------
% Decide if using or not the RISE print and plot settings
USE_RISE_PLOT = true;
USE_RISE_PRINT = false;

%--------------------------------------------------------------------------

%  %--------------------------------------------------------------------------
%  % format of numbers on MATLAB terminal
%  format long g

rise_data=cell(0,2);


if USE_RISE_PRINT
    
    % ------------------------------------------------------------------------
    % Version Variables
    
    NOT_INSTALLED = 'Not installed';
    matlab_version = NOT_INSTALLED;
    symbolic_math_version = NOT_INSTALLED;
    optimization_version = NOT_INSTALLED;
    statistics_version = NOT_INSTALLED;
    rise_version = NOT_INSTALLED;
    
    vs = ver;
    for jj = 1:length(vs)
        v = vs(jj);
        switch v.Name
            case 'MATLAB'
                matlab_version = [v.Version ' ' v.Release];
            case 'Symbolic Math Toolbox'
                symbolic_math_version = [v.Version ' ' v.Release];
            case 'Optimization Toolbox'
                optimization_version = [v.Version ' ' v.Release];
            case 'Statistics Toolbox'
                statistics_version = [v.Version ' ' v.Release];
            case 'RISE Toolbox'
                rise_version = [v.Version ' ' v.Release];
        end
    end
    rise_data=[rise_data
        'matlab_version', matlab_version
        'symbolic_math_version', symbolic_math_version
        'optimization_version', optimization_version
        'statistics_version', statistics_version
        'rise_version', rise_version
        'rise_required_matlab_version', '7.11'
        ];
    
    
    %%  %--------------------------------------------------------------------------
    %%  % Check and load user parameters
    %%  %
    %%  loadPrefs;
    
    %--------------------------------------------------------------------------
    % set page properties for printing
    set(0, 'DefaultFigurePaperOrientation','landscape');
    set(0, 'DefaultFigurePaperType','A4');
    set(0, 'DefaultFigurePaperUnits', 'centimeters');
    set(0, 'DefaultFigurePaperPositionMode', 'manual');
    set(0, 'DefaultFigurePaperPosition', [3.56 2.03 22.56 16.92]);
end


%--------------------------------------------------------------------------
% Plot settings
if USE_RISE_PLOT
    
    % ------------------------------------------------------------------------
    % General Variables
    
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
        'rise_default_plot_colors',{rise_default_plot_colors}];
    
    
    %     set(0, 'DefaultAxesXColor', [0 0 0]);
    %     set(0, 'DefaultAxesYColor', [0 0 0]);
    %    set(0, 'defaultfigurenumbertitle', 'on');
    set(0, 'DefaultFigureColor', 'w');
    %     set(0, 'DefaultFigurePosition', [0 0 1200 700]);
    %     set(0, 'DefaultAxesPosition', [0.13 0.15 0.775 0.75]);
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
rise_root=strrep(which('rise'), fullfile('rise', 'classes', '@rise', 'rise.m'), '');
rise_data=[rise_data
    {'rise_root',rise_root}];
rise_data=[rise_data
    {'rise_pdflatex',~isempty(pdflatex_path)}];

for id=1:size(rise_data,1)
    if flag
        rmappdata(0,rise_data{id,1});
    else
        % Set RISE Root dir
        setappdata(0, rise_data{id,1}, rise_data{id,2});
    end
end

if flag
    set(0,'default');
end
% Show logo
welcome_message();

    function welcome_message()
        
        vv = ver('RISE_Toolbox');
        l1 = '+----------------------------------------------+';
        
        disp(l1);
        disp(['Welcome to the ', vv.Name])
        disp(['Version: ', vv.Version])
        disp(['Release: ', vv.Release])
        disp(['Date: ', vv.Date])
        disp('Please address concerns, problems, suggestions ')
        disp('and desideratas to junior.maih@gmail.com')
        disp('Thank you for your feedback !!!')
        if retcode
            disp('pdflatex/Miktex could not be located')
        end
        disp(l1);
        
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

