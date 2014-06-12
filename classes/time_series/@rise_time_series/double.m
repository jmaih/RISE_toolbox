function out=double(this)
if isempty(this)
    out=[];
else
    out=cell2mat(this.data(2:end,2:end,:));
end
% cell2mat does not render the correct size when input is
% empty. So I need to protect against that somehow as I did in
% the past. But I will take this issue later, maybe
end
