function this=subsref(obj,s)
% H1 line
%
% ::
%
%
% Args:
%
% Returns:
%    :
%
% Note:
%
% Example:
%
%    See also:

% db=ts('1990q1',randn(10,4,3),{'v1','v2','v3','v4'});
% pick a variable:
%                db{'v1'}, db('v1')
% pick variables:
%                db{'v1,v4'}, db('v1,v4') db{{'v1','v4'}}, db({'v1','v4'})
%                db(:,1),
% pick a year :
%                db{'1991'} db('1991')
% pick a date :
%                db{'1992Q1'} db('1992Q1')
% pick dates :
%                db{'1990Q1:1990Q3,1991,1992Q2'}
%                db('1990Q1:1990Q3,1991,1992Q2')
%                db{date2serial('2011Q1'):date2serial('2014Q1')}
%                db(date2serial('2011Q1'):date2serial('2014Q1'))
% pick date and variables :
%                db{'1991Q1','v3'} db('1991Q1','v3')
%                db{'1991Q1',{'v1','v2'}} db('1991Q1',{'v1','v2'})
%                db{'1991:1992',{'v1','v2'}} db('1991:1992',{'v1','v2'})
%                db{'1991Q1:1992Q2',{'v1','v2'}} db('1991Q1:1992Q2',{'v1','v2'})
%                db{:,{'v1','v2'}} db(:,{'v1','v2'})
%                db{'1991Q1:1992Q2',:} db('1991Q1:1992Q2',:)
%                db{:,:} db(:,:)
%                db{date2serial('2011Q1'):date2serial('2014Q1'),{'v1','v2'}}
% pick date and variables and pages:
%                db{'1991Q1',{'v1','v2'},3:4} db('1991Q1',{'v1','v2'},3:4)
%                db{'1991Q1',{'v1','v2'},:} db('1991Q1',{'v1','v2'},:)
%                db{:,:,:} db(:,:,:)
% shifting/lagging: lag or lead series
%                db{-1} db(-1)
%                db{+3} db(+3)
% boolean indexing
%                db{db{'v1'}>0} db{db{'v1'}>0}
%--------------------------------------------------------------------------
% truncate : same as window ?
% plot(db,'subplots',true,'figsize',[3,3],'nticks',7,'logy',true)
% plot(db,'secondary_y',{'V1','V2'})
% bar(db)
% bar([db1,db2])
% barh(db)
% hist(db)
% scatter_matrix([db1,db2],'figsize',[3,3],'diagonal','kde')
% lag_plot : plot(y_t,y_{t+k})
% autocorrelation_plot : plot(y_t): p. 437
% head(db,n) n=min(5,smpl) by default : use date numbers
% tail(db,n) n=min(5,smpl) by default
% index(db) to display the basic infos
% describe(db): cell array: count, mean, std, min,
%                25%, 50%, 75%, max
% transpose(db)

switch s(1).type
    
    case '.'
        % subs = character string
        this=builtin(mfilename,obj,s);
        
    case {'()','{}'}
        [date_numbers,datta,...
            rows_dates,varloc,pages]=process_subs(obj,s(1).subs,mfilename);
        
        if isvector(rows_dates)
            
            rows_dates=rows_dates(:);
            
        end
        
        if size(rows_dates,2)>1||size(rows_dates,3)>1
            
            error('I don''t know how to concatenate multiple series with different dates')
        
        end
        
        date_numbers=date_numbers(rows_dates);
        
        if any(isnan(varloc))
            
            datta=datta(rows_dates,:,pages);
            
            this=ts(date_numbers(:),datta,obj.varnames,obj.description);
            
        elseif isempty(varloc)
            
            this=ts.empty(0);
            
        else
            
            datta=datta(rows_dates,varloc,pages);
            
            this=ts(date_numbers(:),datta,obj.varnames(varloc),...
                obj.description(varloc));
            
        end
        
        s=s(2:end);
        
        if ~isempty(s)
            
            this=subsref(this,s);
            
        end
        
    otherwise
        
        error(['unexpected type "',s(1).type,'"'])
        
end


end
