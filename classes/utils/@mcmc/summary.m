function [myMeanStdev,myQuantiles]=summary(obj,varargin)

defaults={
    'percnt',[2.5,25,50,75,97.5],@(x)all(isnumeric(x))&& all(x>0) && all(x<100)
    'batch_size',100,@(x)isnumeric(x) && isscalar(x) && x>=1 && floor(x)==ceil(x)
    };

[percnt,batch_size]=parse_arguments(defaults,varargin{:});

mycell=cell(obj.nparams,5);

mycell(:,1)=obj.pnames(:);

npc=numel(percnt);

mypercentiles=cell(obj.nparams,npc+3);

mypercentiles(:,1)=obj.pnames(:);

for ipar=1:obj.nparams
    
    x=load_draws(obj,obj.pnames{ipar},[]);
    
    xx=x(:,1:batch_size:end);
    
    xx=xx(:);
    
    x=x(:);
    
    theMean=mean(x);
    
    theStd=std(x);
    
    theStd_batch=std(xx);
    
    mycell{ipar,2}=theMean;
    
    mycell{ipar,3}=theStd;
    
    mycell{ipar,4}=theStd_batch;
    
    mycell{ipar,5}=theMean/theStd;
    
    mypercentiles{ipar,2}=obj.best.x(ipar);
    
    mypercentiles{ipar,3}=median(x);
    
    mypercentiles(ipar,4:end)=num2cell(prctile(x,percnt));
    
end

mycell=[{'','mean','SD','SD(batch)','mean/SD'}
    mycell];

num2cell_percnt=num2cell(percnt);

for ii=1:numel(percnt)
    
    num2cell_percnt{ii}=sprintf('%0.1f%%',num2cell_percnt{ii});
    
end

mypercentiles=[{'','mode','median'},num2cell_percnt
    mypercentiles];

if nargout
    
    myMeanStdev=mycell;
    
    myQuantiles=mypercentiles;
    
else
    
fprintf(1,'Iterations = 1:%0.0f\n',obj.npop);

fprintf(1,'Thinning interval = %0.0f\n',nan);

fprintf(1,'Number of chains = %0.0f\n',obj.nchains);

fprintf(1,'Sample size per chain = %0.0f\n',obj.npop);

fprintf(1,'Percentage of draws discarded/dropped = %0.4f percent\n',100*obj.i_dropped);

fprintf(1,'\nEmpirical means and standard deviations\n');

disp(mycell)

fprintf(1,'Quantiles for each parameter\n');

disp(mypercentiles)

end

end
