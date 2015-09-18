function [R,HHRestr]=build_shock_restrictions(H,G,yrest_id,xrest_id,...
    nsteps)

[model,ny]=build_model();
y0=zeros(ny,1);

[M,~,~,~,HHRestr]=utils.forecast.rscond.stochastic_impact_cumulator(model,...
    y0,nsteps,yrest_id,xrest_id);

R=[M.R;M.S];

    function m=build_model()
        ny=size(H,1);
        T=[H,zeros(ny,1),G(:,:)];
        sstate=zeros(ny,1);
        state_cols=1:ny;
        k=size(G,3)-1;
        Qfunc=@(x)1;
        m=struct('T',{T},'sstate',{sstate},'state_cols',state_cols,...
            'k',k,'Qfunc',Qfunc);
    end
end