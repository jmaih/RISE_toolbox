function [yh,res]=chowlin(y0,Xh0,varargin)
% CHOWLIN - Temporal disaggregation using the Chow-Lin method
%
% ::
%
%
%   [yh,res]=chowlin(y0,Xh0)
%   [yh,res]=chowlin(y0,Xh0,aggreg_type)
%   [yh,res]=chowlin(y0,Xh0,aggreg_type,estim_method)
%   [yh,res]=chowlin(y0,Xh0,aggreg_type,estim_method,ngrid)
%
% Args:
%
%    - **y0** [ts] : low-frequency left-hand-side variable
%      the time series. The following are admitted:
%      - annual data : e.g. 1990 or '1990'
%      - bi-annual data : e.g. '1990H1'
%      - Quarterly data : e.g. '1990Q3'
%      - monthly data : e.g. '1990M12'
%    - **Xh0** [ts|struct] : high-frequency right-hand-side explanatory
%          variables
%    - **aggreg_type** ['flow'|{'average'}|'index'|'last'|'first'] : type of
%      aggregation:
%      - 'flow' (or 1) is the sum,
%      - 'average','index' (or 2) is the average,
%      - 'last' (or 3) is the last element,
%      - 'first' (or 4) is the first element,
%    - **estim_method** [{0}|1|2] : estimation method
%      - 0 : Generalized/Weighted Least squares (with grid)
%      - 1 : Maximum Likelihood (with grid),
%      - 2 : Maximum Likelihood (without grid) using fmincon as optimizer
%    - **ngrid** [integer|{250}] : number of grid points
%
% Returns:
%    :
%
%    - **yh** [ts] : (disaggregated) high-frequency left-hand-side variable
%
%    - **res** [struct] : structure with further details on the computations
%
% Note:
%
% Example:
%
%    See also:
%    ---------
%
%    REFERENCES:
%    -----------
%    - Chow, G. and Lin, A.L. (1971) "Best linear unbiased
%    distribution and extrapolation of economic time series by related
%    series", Review of Economic and Statistics, vol. 53, n. 4, p. 372-375.
%    - Bournay, J. y Laroque, G. (1979) "Reflexions sur la methode
%    d'elaboration  des comptes trimestriels", Annales de l'INSEE, n. 36, p.
%    3-30.

check_inputs()

dnl=y0.date_numbers;

dnh=Xh0.date_numbers;

fh=utils.time_series.freq2freq(get(Xh0,'frequency'));

fl=utils.time_series.freq2freq(get(y0,'frequency'));

s=fh/fl;

y=double(y0);

[startl,starth]=find_start();

if startl~=1
    
    dnl=y0.date_numbers;
    
    start_date=parser.any2str(serial2date(dnl(startl)));
    
    error(['Start date for low frequency does not match the high-frequency ',...
        'information. You may want to discard the first ',int2str(startl-1),...
        ' observations in the low-frequency database or more specifically ',...
        'start the low-frequency database in ',start_date])
    
end

Xh=double(Xh0);

% Low frequencies without corresponding high frequency discarded
%----------------------------------------------------------------
Xh(1:starth-1,:)=[];

if any(isnan(y0(:)))||any(isnan(Xh(:)))
    
    error('nan in the data')
    
end

res=utils.time_series.mychowlin(y,Xh,s,varargin{:});

yh=reset_data(Xh0,res.y,'');

    function [startl,starth]=find_start()
        
        [decl,~]=decompose_date(y0.date_numbers);
        
        [dech,~]=decompose_date(Xh0.date_numbers);
        
        year_l=[decl.year];
        
        Period_l=[decl.period]*s;
        
        year_h=[dech.year];
        
        Period_h=[dech.period];
        
        if (year_l(end)>year_h(end))||...
                (year_l(end)==year_h(end) && ...
                Period_l(end)>Period_h(end))
            
            end_datel=parser.any2str(serial2date(dnl(end)));
    
            end_dateh=parser.any2str(serial2date(dnh(end)));
    
            error(['Low-frequency database cannot end (',end_datel,')',...
                ' after the high-frequency database (',end_dateh,')'])
            
        end
            
        
        feasible=false;
        
        startl=0;
        
        while ~feasible
            
            startl=startl+1;
            
            ystart=find(year_h==year_l(startl) & Period_h==Period_l(startl),1);
            
            if ~isempty(ystart)
                
                starth=find(year_h==year_h(ystart) & ...
                    Period_h==Period_h(ystart)-s+1, 1);
                
                feasible=~isempty(starth);
            
            end
            
        end
        
    end

    function check_inputs()
        
        if ~isa(y0,'ts')||get(y0,'NumberOfVariables')~=1
            
            error('first input must be a scalar ts object')
            
        end
        
        if ~isa(Xh0,'ts')
            
            Xh0=ts.collect(Xh0);
            
        end
        
        n=length(varargin);
        
        if n==0,return,end
        
        if ~isempty(varargin{1})
            
            aggreg_type={'flow','average','index','last','first',1,2,3,4};
            
            v3=varargin{1};
            
            cond1=ischar(v3) && ismember(v3,aggreg_type(1:5));
            
            cond2=isnumeric(v3) && any(v3==cell2mat(aggreg_type(6:end)));
            
            if ~(cond1||cond2)
                
                disp(aggreg_type)
                
                error('Third input (aggreg_type) must be one of the above')
                
            end
            
        end
        
        if n==1,return,end
        
        if ~isempty(varargin{2})
            
            estim_method=0:2;
            
            if ~ismember(varargin{2},estim_method)
                
                error('Fourth input (estim_method) must be 0, 1 or 2')
                
            end
            
        end
        
        if n==2,return,end
        
        if ~isempty(varargin{3})% ngrid
            
            v5=varargin{3};
            
            if ~(isa(v5,'double') && v5>0 && isfinite(v5) && ceil(v5)==floor(v5))
                
                error('Fifth input (ngrid) must be a finite and positive integer')
                
            end
            
        end
        
        if n>3
            
            error('too many input arguments')
            
        end
        
    end

end