function keep(varargin)
% keep - clears the workspace for undesired variables
%
% ::
%
%
%   keep v1 v2 v3...
%
%   keep('v1','v2','v3',...,'vn')
%
% Args:
%
%    - **v1, v2,..,vn** [char]: names of the variables to clear
%
% Returns:
%    :
%
%    - none
%
% Note:
%
% Example:
%
%    See also: clear

% This is a modification of keep.m by Xiaoning (David) Yang xyang@lanl.gov
% 1998.

% Keep all
%----------
if isempty(varargin)
    return
end

% See what are in caller workspace
%---------------------------------
wh = evalin('caller','who');

% Check workspace variable
%--------------------------
if isempty(wh)
    return
end

myvars=varargin;
if iscell(myvars{1})
    if numel(myvars)>1
        error('inputs must char or a cell array')
    end
    myvars=myvars{1};
end

% check that the variables to keep are all present
%--------------------------------------------------

n=numel(myvars);
present=false(1,n);
for iv=1:n
    present(iv)=any(strcmp(myvars{iv},wh));
end

if any(~present)
    disp(myvars(~present))
    error('The variables above are not present in the workspace')
end

discard=setdiff(wh,myvars);

mystring=cell2mat(strcat('''',discard(:)',''','));

evalin('caller',['clear(',mystring(1:end-1),')'])
