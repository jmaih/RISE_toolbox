function new_date=obs2date(start_date,obs)
start_serial=date2serial(start_date);
newobs=1;
while newobs<obs
    newobs=newobs+1;
end
new_serial=start_serial+newobs-1;
new_date=serial2date(new_serial);
end

