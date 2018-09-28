
nreplic=1000;

N=100000;

order=6;

for k=1:5 % degrees of freedom of the chi-square
    
    % theoretical moments
    theoMoms=moments.tailored_chisquare(order,k);
    
    empMoms=zeros(nreplic,order);
    
    for irep=1:nreplic
        
        % draws
        d0=randn(N,k); d1=sum(d0.^2,2); d2=(d1-k)/sqrt(2*k);
        
        % empirical moments
        empMoms(irep,:)=sum(d2.^(1:order),1)/N;
        
    end
    
    figure('name',['comparing empirical vs theoretical moments for k=',int2str(k)]);
    
    for io=1:order
        
        [F,XI]=ksdensity(empMoms(:,io));
        subplot(3,2,io)
        plot(XI,F,'linewidth',2)
        hold on
        stud=F==max(F);
        plot(XI(stud)*ones(1,2),[0,F(stud)],'linewidth',2)
        axis tight
        title({['m',int2str(io)],['theo=',num2str(theoMoms(io),2),...
            ' -- emp=',num2str(XI(stud),2)]})
        
    end
    
end

%%

z=zeros(N+1,1);

rho=0.8;

sigma=0.05;

for ii=2:N+1
   
    z(ii)=rho*z(ii-1)+sigma*d2(ii-1);
    
end

close all

figure()
subplot(2,1,1)
plot(d2),hold on, plot([1,N],[0,0],'linewidth',2)
axis tight

subplot(2,1,2)
plot(z),hold on, plot([1,N],[0,0],'linewidth',2)
axis tight