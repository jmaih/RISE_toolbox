function [nticks,varargs]=extract_number_of_ticks(varargin)

nticks=[];
n=length(varargin);
extract=[];
for icol=1:n
    if strcmpi(varargin{icol},'nticks')
        nticks=varargin{icol+1};
        extract=[icol,icol+1];
        break
    end
end
if ~isempty(extract)
    varargin(extract)=[];
end
varargs=varargin;