function out=if_elseiff(args,varargin)
% every second element of varargin is either:
% - a function handle, in which case it is evaluated as a function
% - anything else, in which case the entry is returned as is
%
% abs(if_elseiff({1},true,@cos,false,@sin)-0.54030230586814)<1e-15
% abs(if_elseiff({1},false,@cos,true,@sin)-0.841470984807897)<1e-15
%
% abs(if_elseiff({1},true,3.5,false,7.2)-3.5)==0
% abs(if_elseiff({1},false,3.5,true,7.2)-7.2)==0

stud=cell2mat(varargin(1:2:end));

funcs=varargin(2:2:end);

ping=find(stud,1,'first');

out=funcs{ping};

if isa(out,'function_handle')
    
    out=out(args{:});
    
end

end

