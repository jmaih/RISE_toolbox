function [Hist_db,Cond_db]=xl2databases(xlsfilename,varobs,condvarobs,horizon,startdate,enddate,sheet,range)
%{
startdate='1947q3';
enddate='2004Q4';
varobs=char('dy','dc','dinve','labobs','pinfobs','dw','robs');
condvarobs=char('dy','dc','pinfobs');
[Hist_db,Cond_db]=Historical_and_Conditional_databases('usmodel_data_spf',varobs,condvarobs,6,startdate,enddate)
%}

disp([mfilename,':: ',upper('this file may need updating... contact junior.maih@gmail.com if it does not work')])

if nargin<8
    range=[];
    if nargin<7
        sheet=[];
        if nargin<6
            enddate='';
            if nargin<5
                startdate=nan;
                if nargin<4
                    horizon=0;
                    if nargin<3
                        condvarobs='';
                        if nargin<2
                            error([mfilename,':: At least the xls file and the list of variables should be provided'])
                        end
                    end
                end
            end
        end
    end
end
if ~isempty(range)
    rangeflag=true;
end
if ~isempty(sheet)
    sheetflag=true;
end
if sheetflag && rangeflag
    [datta,txt]=xlsread(xlsfilename,sheet,range);
elseif rangeflag && ~sheetflag
    error([mfilename,':: Range cannot be provided without the sheet'])
elseif ~rangeflag && sheetflag
    [datta,txt]=xlsread(xlsfilename,sheet);
else
    [datta,txt]=xlsread(xlsfilename);
end


VarNames=char(txt(1,2:end));
VarNames_id=[];
for j=1:size(VarNames,1)
    vj=deblank(VarNames(j,:));
    vj_id=strmatch(vj,VarNames,'exact');
    if isempty(vj)||numel(vj_id)>1
        continue
    end
    VarNames_id=[VarNames_id,vj_id]; %#ok<AGROW>
end
VarNames=VarNames(VarNames_id,:);

% locate the column with the dates and get the date
if ~isempty(txt{2,1})
    begin=txt{2,1};
    datta=datta(:,VarNames_id);
else% then the dates are the first column of the data
    begin=num2str(datta(1,1));
    begin=[begin(1:4),'Q',begin(end)];
    datta=datta(:,VarNames_id+1);
end

if all(isnan(startdate))||isempty(startdate)
    startdate=begin;
end

% create a raw database
DB0=ts(begin,datta,VarNames);

% restrict attention to the relevant sample
DB0=DB0.window(startdate,enddate);

% extract observer equation variables
Hist_db = DB0.window(startdate,enddate,varobs);

% build conditional database
datta=cell2mat(DB0.data(2:end,2:end));
nvars=size(condvarobs,1);
Varscond=nan(nvars,horizon);
for v=1:nvars
    vj=deblank(condvarobs(v,:));
    for h=1:horizon
        vjh_id=strmatch([vj,int2str(h)],DB0.varnames,'exact');
        if ~isempty(vjh_id)
            Varscond(v,h)=vjh_id;
        else
            disp got one
        end
    end
end

NumberOfObservations=size(datta,1);
datamat=nan(horizon,nvars,NumberOfObservations);
for o=1:NumberOfObservations
    for v=1:nvars
        good=~isnan(Varscond(v,:));
        datamat(good,v,o)=datta(o,Varscond(v,good))';
    end
end
Cond_db=ts(startdate,permute(datamat,[3,2,1]),condvarobs);