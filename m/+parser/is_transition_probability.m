function [istp,diagonal,chain_name,max_state]=is_transition_probability(name)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% Description
% ------------
%
% Examples
% ---------
%
% See also: 

istp=false;
diagonal=false;
chain_name='';
max_state=[];
if isletter(name(1))
    underscore=strfind(name,'_');
    if ~isempty(underscore)&& numel(underscore)==3
        if size(name,2)>=8 && strcmp('tp',name(underscore(1)+1:underscore(2)-1))
            first=name(underscore(2)+1:underscore(3)-1);
            second=name(underscore(3)+1:end);
            if ~all(isletter(first)) && ~all(isletter(second))
                first=str2double(first);
                second=str2double(second);
                if ~isnan(first) && ~isnan(second)
                    istp=true;
                    chain_name=name(1:underscore(1)-1);
                    max_state=max(first,second);
                    if second==first
                        diagonal=true;
                    end
                end
            end
        end
    end
end