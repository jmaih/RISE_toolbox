function rise_startup(flag)
 if nargin<1
	 flag=false;
 end
 testing=true;

 % this function sets up RISE (it replaces setpaths)
  
  %-----------------------------------------------------------------------
  % Decide if using or not the RISE print and plot settings
  USE_RISE_PLOT = false;
  USE_RISE_PRINT = true;
  
  %--------------------------------------------------------------------------  
  
  %--------------------------------------------------------------------------
  %--------------------------------------------------------------------------
  % NO NEED TO EDIT BELOW HERE
  %--------------------------------------------------------------------------
  %--------------------------------------------------------------------------
  
  v = ver('RISE');
  
%  %--------------------------------------------------------------------------
%  % format of numbers on MATLAB terminal
%  format long g
  

  
  % ------------------------------------------------------------------------
  % General Variables
  setappdata(0, 'xmlsetsize', 50000); % Max size of an xml data set <Set></Set>
  
  setappdata(0, 'rise_default_plot_colors', { ...
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
    [0.9 0.266 0.593]});
  
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
  
  setappdata(0, 'matlab_version', matlab_version);
  setappdata(0, 'symbolic_math_version', symbolic_math_version);
  setappdata(0, 'optimization_version', optimization_version);
  setappdata(0, 'statistics_version', statistics_version);
  setappdata(0, 'rise_version', rise_version);
  setappdata(0, 'rise_required_matlab_version', '7.11');
  
  
%%  %--------------------------------------------------------------------------
%%  % Check and load user parameters
%%  %
%%  loadPrefs;
    
  %--------------------------------------------------------------------------
  % set page properties for printing
  if USE_RISE_PRINT
    set(0, 'DefaultFigurePaperOrientation','landscape');
    set(0, 'DefaultFigurePaperType','A4');
    set(0, 'DefaultFigurePaperUnits', 'centimeters');
    set(0, 'DefaultFigurePaperPositionMode', 'manual');
    set(0, 'DefaultFigurePaperPosition', [3.56 2.03 22.56 16.92]);
  end
  
%%%  % ------------------------------------------------------------------------
%%%  % Backup MATLAB's plot settings
%%%  utils.plottools.backupDefaultPlotSettings();
  
  %--------------------------------------------------------------------------
  % Plot settings
  if USE_RISE_PLOT
    set(0, 'DefaultAxesXColor', [0 0 0]);
    set(0, 'DefaultAxesYColor', [0 0 0]);
%    set(0, 'defaultfigurenumbertitle', 'on');
    set(0, 'DefaultFigureColor', 'w');
    set(0, 'DefaultFigurePosition', [0 0 1200 700]);
    set(0, 'DefaultAxesPosition', [0.13 0.15 0.775 0.75]);
  end
  
  % Add user model paths
%%%  prefs = getappdata(0, 'RISEpreferences');
%%%  searchPaths = prefs.getModelsPrefs.getSearchPaths;
  searchPaths=collect_paths(mfilename);
  for jj = 1:numel(searchPaths)
  	if flag
	    rmpath(searchPaths{jj});
	else
	    addpath(searchPaths{jj});
	end
  end
  
%  % Install extensions
%  utils.helper.installExtensions;
  
  %--------------------------------------------------------------------------
  % Activate correct helptoc.xml file (depending on MATLAB version)

  if ~testing
	  % Define MATLAB helptoc version
	  matlabRelease = version('-release');
	  switch matlabRelease
	    case {'2008a', '2008b', '2009a'}
	      matlabRelease = 'R2009a';
	    case '2009b'
	      matlabRelease = 'R2009b';
	    case '2010a'
	      matlabRelease = 'R2010a';
	    otherwise
	      matlabRelease = 'R2010a';
	  end
	  
	  % Get info.xml path from the path of rise_startup because it might be
	  % happen that the RISE toolbox is not on top of the MATLAB path
	  startupPath = fileparts(which(mfilename()));
	  infoPath = strrep(startupPath, strcat('m', filesep(), 'etc'), '');
	  
	  % read info.xml file in order to get the helptoc.xml path
	  infoXML = xmlread(fullfile(infoPath, 'info.xml'));
	  tbNameNode = infoXML.getElementsByTagName('name');
	  tbName = tbNameNode.item(0).getFirstChild.getData;
	  if strcmp(tbName, 'RISE')
	    helpLocationNodes = infoXML.getElementsByTagName('help_location');
	    helpLocation = char(helpLocationNodes.item(0).getFirstChild.getTextContent);
	  else % Otherwise error out
	    error('Can not find info.xml file for My Toolbox');
	  end
	  
	  helptocLocation = fullfile(infoPath, helpLocation);
	  
	  helptocSource = fullfile(helptocLocation, strcat('helptoc', matlabRelease, '.xml'));
	  helptocDest   = fullfile(helptocLocation, 'helptoc.xml');
	  
	  copyfile(helptocSource, helptocDest);
  end
  
  % Set RISE Root dir
  riseroot = strrep(which('rise'), fullfile('rise', 'classes', '@rise', 'rise.m'), '');
  setappdata(0, 'RISEROOT', riseroot);
  
  
  % Show logo
  showLogo();
      
end


function showLogo()
  
  v = ver('RISE_Toolbox');
  
  logo = {...
    '                                        ',...
    '                  ****                  ',...
    '                   **                   ',...
	'  _____	   _    _____ 	 _____ 		 ',...
	' |  _  \	  |*|  | ____|  | ____|   	 ',...
	' | |_| |	   _   | |      | |       	 ',...
	' |     /	  |	|  | |__ _  | |___    	 ',...
	' | |\ \	  |	|  |____  | |  ___|   	 ',...
	' | | \ \	  | |       | | | |       	 ',...
	' | |  \ \    | |   ____| | | |___       ',...
	' |_|	\_\   |_|  |______| |_____|   	 ',...
    '                   **                   ',...
    '                  ****                  ',...
    };
  
  l1 = '+----------------------------------------------------+';
%   ll = length(l1);
  
  disp(l1);
%   for jj = 1:numel(logo)
%     fprintf(1,'%s\n',logo{jj});
%   end
  disp(['Welcome to the ', v.Name])
  disp(['Version: ', v.Version])
  disp(['Release: ', v.Release])
  disp(['Date: ', v.Date])
  disp(l1);
  
end


function list=collect_paths(filename)
fullpath=which(filename);
loc=strfind(fullpath,filename);
tmp=genpath(fullpath(1:loc-2));

if ispc
    semicols=strfind(tmp,';');
elseif ismac
    semicols=strfind(tmp,':');
else
    error([mfilename,':: unknown system '])
end
previous=0;
collect={};
for ii=1:numel(semicols)
    current=semicols(ii);
    thepath=tmp(previous+1:current-1);
    if isempty(strfind(thepath,'.svn')) && ...
            isempty(strfind(lower(thepath),'test')) && ...
            isempty(strfind(lower(thepath),'junk'))
       collect=[collect,{thepath}];
    end
    previous=current;
end
if nargout
	list=collect;
end

end

% END
