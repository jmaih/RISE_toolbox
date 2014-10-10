function tmp=main_frame(self,name_down)
% H1 line
%
% Syntax
% -------
% ::
%
% Inputs
% -------
%
% Outputs
% --------
%
% Description
% ------------
%
% Examples
% ---------
%
% See also: 

if nargin<2
    name_down=true;
end

if ~isempty(self.date_numbers)
    if ~all(cellfun(@isempty,self.varnames))
        Header=[{'Time'};self.varnames(:)];
    else
        Header=[{'Time'};repmat({''},self.NumberOfVariables,1)];
    end
    Header=Header';
    data=self.data;
    [these_dates,this.frequency]=serial2date(self.date_numbers);
    if self.cell_style
        nrows=self.NumberOfVariables+1;
        tmp=cell(nrows,nrows);
        tmp(1,:)=Header;
        tmp(:,1)=Header';
        tmp=repmat(tmp,[1,1,self.NumberOfPages]);
        these_dates=repmat(these_dates,[1,1,self.NumberOfPages]);
        tmp0=tmp;
        tmp=repmat(tmp0,[self.NumberOfObservations,1,1]);
        offset=0;
        for iobs=1:self.NumberOfObservations
            tmp0(1,1,:)=these_dates(iobs,1,:);
            tmp0(2:end,2:end,:)=num2cell(data{iobs});
            tmp(offset+(1:nrows),:,:)=tmp0;
            offset=offset+nrows;
        end
    else
        tmp=cell(self.NumberOfObservations+1+name_down,self.NumberOfVariables+1);
        tmp(1,:)=Header;
        if name_down
            tmp(end,:)=Header;
        end
        tmp(2:end-name_down,1)=these_dates;
        tmp=repmat(tmp,[1,1,self.NumberOfPages]);
        if ~isempty(data)
            for ipage=1:self.NumberOfPages
                tmp(2:end-name_down,2:end,ipage)=num2cell(data(:,:,ipage));
            end
        end
    end
end
