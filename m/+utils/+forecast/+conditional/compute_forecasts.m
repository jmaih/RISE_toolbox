function Y=compute_forecasts(Y,H,G,E,Nsteps,NumberOfAnticipatedPeriods)
if nargin<6
    NumberOfAnticipatedPeriods=0;
end
% Compute conditional forecasts
for t=1:max(0,Nsteps)
    Y(:,t+1)=H*Y(:,t);
    for j=1:NumberOfAnticipatedPeriods
        Y(:,t+1)=Y(:,t+1)+G(:,:,j)*E(:,t+j-1);
    end
end
end

