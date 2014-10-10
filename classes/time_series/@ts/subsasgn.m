function obj=subsasgn(obj,s,b)% subsasgn
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


switch s.type
    case '.'
        % subs = character string
        obj=builtin(mfilename,obj,s,b);
    case {'()','{}'}
        [date_numbers,datta,...
            rows_dates,varloc,pages]=process_subs(obj,s.subs,mfilename);
        %         date_numbers=date_numbers(rows_dates);
        if isvector(rows_dates)
            if isscalar(b)
                datta(rows_dates,varloc,pages)=b;
            else
                datta(rows_dates,varloc,pages)=b;
            end
        else
            datta(rows_dates)=b;
        end
%         datta=datta(rows_dates,varloc,pages);
        obj=ts(date_numbers(:),datta,obj.varnames);
    otherwise
        error(['unexpected type "',s.type,'"'])
end

end
