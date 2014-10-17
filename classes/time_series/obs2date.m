function new_date=obs2date(start_date,obs)
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
% More About
% ------------
%
% Examples
% ---------
%
% See also: 

start_serial=date2serial(start_date);
nobs=numel(obs);
new_serial=nan(size(obs));
pos=obs>0;
neg=obs<0;
if any(pos)
    new_serial(pos)=do_positives(obs(pos));
end
if any(neg)
    new_serial(neg)=do_negatives(obs(neg));
end
if any(~pos & ~neg)
    error('observation numbers can only be positive on negative integers')
end

new_date=serial2date(new_serial);

    function new_serial=do_positives(obs)
        newobs=1*ones(size(obs));
        for iobs=1:nobs
            while newobs(iobs)<obs(iobs)
                newobs(iobs)=newobs(iobs)+1;
            end
        end
        new_serial=start_serial+newobs-1;
    end

    function new_serial=do_negatives(obs)
        new_serial=start_serial+obs;
    end
end

