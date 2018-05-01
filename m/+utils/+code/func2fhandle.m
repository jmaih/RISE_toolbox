function varargout = func2fhandle(varargin)
% func2fhandle - creates a handle to an m-file not on the matlab search path
%
% ::
%
%
%   hdl = func2fhandle(full_path/func)
%   hdl = func2fhandle(relative_path/func)
%   [hdl1,hdl2,...,hdln] = func2fhandle(p1,p2,...,pn)
%   [hdl1,hdl2,...,hdln] = func2fhandle({p1,p2,...,pn})
%
% Args:
%
%    - **pth** [char|cellstr] : path names to create handles for. The path
%      could be relative to the current directory or absolute
%
% Returns:
%    :
%
%    - **hdl** [function_handle] : function handle(s)
%
% Note:
%
% Example:
%
%    See also:

% restore the current directory when we are done
currentDir = pwd;
c = onCleanup(@()cd(currentDir));

if iscell(varargin{1})
    if numel(varargin)>1
        error('wrong input specification')
    end
    varargin=varargin{1};
end
varargout=varargin;

nf=length(varargin);

for ifunc=1:nf
    fff=varargin{ifunc};
    if ~isa(fff, 'function_handle')
        if ~ischar(fff)
            error('all inputs must be function handles or char')
        end
        if strcmp(fff(end-1:end),'.m')
            fff=fff(1:end-2);
        end
        if ~exist([fff,'.m'],'file')
            error(['file ',fff,'.m not found'])
        end
        [pth,fname] = fileparts(fff);
        cd(pth)
        varargout{ifunc}= str2func(['@',fname]);
        cd(currentDir)
    end
end

end
