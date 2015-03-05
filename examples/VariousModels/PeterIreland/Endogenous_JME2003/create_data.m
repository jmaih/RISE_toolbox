function data=create_data(estim_start_date,estim_end_date)
if nargin<2
    estim_end_date=[];
    if nargin<1
        estim_start_date=[];
    end
end
smpl_start_date='1959Q1';
mydat=load('cimpr.dat');
bigt = size(mydat,1);
if isempty(estim_start_date)
    first_obs=1;
else
    first_obs=date2obs(smpl_start_date,estim_start_date);
end
if isempty(estim_end_date)
    last_obs=bigt;
else
    last_obs=date2obs(smpl_start_date,estim_end_date);
end

ct=mydat(:,1);
it=mydat(:,2);
mt=mydat(:,3);
trend = (1:bigt)';
xxx = [ones(bigt,1),trend];

xxx=xxx(first_obs:last_obs,:);

betac = xxx\ct(first_obs:last_obs);
betai = xxx\it(first_obs:last_obs);
betam = xxx\mt(first_obs:last_obs);

ct = ct - betac(2)*trend;
it = it - betai(2)*trend;
mt = mt - betam(2)*trend;

mydat(:,1:3)=[ct,it,mt];
data=ts(smpl_start_date,mydat,{'LC','LI','LM','LPI','LR'});
data=pages2struct(data);
data.TREND=ts(data.LC.start,trend);

end