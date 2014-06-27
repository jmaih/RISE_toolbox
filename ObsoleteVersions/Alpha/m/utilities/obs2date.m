function new_date=obs2date(start_date,obs)
start_serial=date2serial(start_date);
if obs>0
    newobs=1;
    while newobs<obs
        newobs=newobs+1;
    end
    new_serial=start_serial+newobs-1;
elseif obs<0
    new_serial=start_serial+obs;
else
    error('observation numbers can only be positive on negative integers')
end
new_date=serial2date(new_serial);
end

