function [look_behind,loc_,look_forward]=look_around(tokk,to_process,preserve_space)
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

%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

if nargin<3
    preserve_space=false;
end
loc_=strfind(to_process,tokk);
look_forward=to_process;
look_behind='';
if ~isempty(loc_)
    loc_=loc_(1);
    look_behind=to_process(1:loc_-1);
    if ~preserve_space
        look_behind(isspace(look_behind))=[];
    end
    span=length(tokk);
    look_forward=to_process(loc_+span:end);
end
end


