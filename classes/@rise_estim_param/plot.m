function plot(obj)
nobj=numel(obj);
N=600;
cutoff=1e-4;
for ii=1:nobj
    ab_space=obj.ab_space;
    figure('name',['pdf + cdf + inv_cdf (',obj(ii).distribution,'(',num2str(obj(ii).a),',',num2str(obj(ii).b),')) for parameter ',obj(ii).name]);
    % functions of the distribution
    [lpdfn,cdfn,icdfn]=distributions.(obj(ii).distribution)();
    subplot(2,2,1)
    xx=linspace(obj(ii).lb,obj(ii).ub,N);
    dens=exp(lpdfn(xx,obj(ii).a,obj(ii).b,obj(ii).c,obj(ii).d));
    significant=dens>cutoff;
    while sum(significant)<N/3
        xx=0.1*xx;
        dens=exp(lpdfn(xx,obj(ii).a,obj(ii).b,obj(ii).c,obj(ii).d));
        significant=dens>cutoff;
    end
    xx=xx(significant);dens=dens(significant);
    plot(xx,dens); axis tight
    title('probability density function (PDF)')
    subplot(2,2,2)
    plot(xx,cdfn(xx,obj(ii).a,obj(ii).b,obj(ii).c,obj(ii).d)); axis tight
    title('cumulative density function (CDF)')
    subplot(2,2,3)
    u=linspace(0,1,N);
    plot(u,icdfn(u,obj(ii).a,obj(ii).b,obj(ii).c,obj(ii).d)); axis tight
    title('Inverse cumulative density function')
    subplot(2,2,4)
    resfunc=@(a,b)sqrt((cdfn(obj(ii).plb,a,b,obj(ii).c,obj(ii).d)-.5*(1-obj(ii).interval_probability)).^2+...
        (cdfn(obj(ii).pub,a,b,obj(ii).c,obj(ii).d)-.5*(1+obj(ii).interval_probability)).^2);
    N=25;
    avals=obj(ii).a;
    bvals=obj(ii).b;
    n=0;
    increment=0.001;
    while n<N
        n=n+1;
        % hyperparameter a
        left=avals(1)-increment;
        right=avals(end)+increment;
        if left<ab_space(1,1)
            left=[];
        end
        if right>ab_space(1,2)
            right=[];
        end
        avals=[left,avals,right];
        % hyperparameter b
        left=bvals(1)-increment;
        right=bvals(end)+increment;
        if left<ab_space(2,1)
            left=[];
        end
        if right>ab_space(2,2)
            right=[];
        end
        bvals=[left,bvals,right];
    end
    [aa,bb]=meshgrid(avals,bvals);
    zz=resfunc(aa,bb);
    surf(aa,bb,zz)
    title('Hyperparameters peak')
    xlabel('a'),ylabel('b'),zlabel('residuals norm')
    supertitle=[obj(ii).name,'::',obj(ii).distribution,'(',num2str(obj(ii).a),',',num2str(obj(ii).b),')'];
    [~,hand]=sup_label(supertitle,'t');set(hand,'FontSize',13)
end
end
