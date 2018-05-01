function varargout=myplot(plotfunc,varargin)
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

funct_type=func2str(plotfunc);

[this,rise_items,matlab_args]=parse_plot_args(varargin{:});
ndatasets=numel(this);
datta=cell(1,ndatasets);
sizdata=cell(1,ndatasets);
date_numbers=cell(1,ndatasets);
for its=1:ndatasets
    datta{its}=this{its}.data;
    date_numbers{its}=this{its}.date_numbers;
    sizdata{its}=size(datta{its});
    if sizdata{its}(2)>1 && numel(sizdata{its})>2 && sizdata{its}(3)>1
        error([mfilename,':: cannot handle many variables and many pages simultaneously'])
    else
        datta{its}=squeeze(datta{its});
    end
    if rise_items.logy
        datta{its}=log(datta{its});
    end
end

[datta,date_numbers,sizdata]=stretch_data(datta,date_numbers,sizdata);

% collect the linewidth info for the horizontal and vertical lines and
% possibly plotyy
%-----------------------------------------------------------------
linewidth_loc=[];
linewidth={};
if ~isempty(matlab_args)
    for iarg=1:numel(matlab_args)
        if ischar(matlab_args{iarg}) && strcmp(matlab_args{iarg},'linewidth')
            linewidth_loc=[iarg,iarg+1];
        end
    end
end
if ~isempty(linewidth_loc)
    linewidth={'linewidth',matlab_args{linewidth_loc(2)}};
    if strcmp(funct_type,'plotyy')
        matlab_args(linewidth_loc)=[];
    end
end
nticks=rise_items.nticks;

% create a unique set of date_numbers and harmonize the datasets
%---------------------------------------------------------------
xdate_numbers=unique([date_numbers{:}]);
if numel(ndatasets)>1
    reset_data(xdate_numbers)
end

xrange=rise_items.xrange;

if ~isempty(xrange)
    reset_data(xrange)
end

if strcmp(funct_type,{'plot_real_time'})
    if ndatasets~=1
        error('number of datasets cannot exceed 1 in real-time plots')
    end
    [nr,nc]=size(datta{1});
    xtime=date_numbers{1}(1):date_numbers{1}(end)+nc+1;
    pp=plot_specs(xtime,nticks,rise_items.date_format);
else
    pp=plot_specs(date_numbers{1},nticks,rise_items.date_format);
end

   [varargout{1:nargout}]=plot_it(datta);
   
   if nargout==0
       clear plot_handle
   end

    function reset_data(xrange)
        for its_=1:ndatasets
            newdata=nan([numel(xrange),sizdata{its_}(2:end)]);
            locs=nan(sizdata{its_}(1),1);
            processed=false(sizdata{its_}(1),1);
            for iloc=1:sizdata{its_}(1)
                ifind=find(date_numbers{its_}(iloc)==xrange);
                processed(iloc)=~isempty(ifind);
                if processed(iloc)
                    locs(iloc)=ifind;
                end
            end
            good=~isnan(locs);
            newdata(locs(good),:,:)=datta{its_}(good,:,:);
            datta{its_}=newdata;
            date_numbers{its_}=xrange;
        end
    end

    function varargout=plot_it(d)
        nout=nargout;
        if any(strcmp(funct_type,{'boxplot'}))
            plotfunc(d{:},this{1}.varnames,matlab_args{:});
            vout={gca()};
        elseif strcmp(funct_type,{'hist'})
            [vout{1:nout}]=plotfunc(d{:},matlab_args{:});
        elseif strcmp(funct_type,{'plot_real_time'})
            ymax=max(d{:}(:));
            ymin=min(d{:}(:));
            map=getappdata(0,'rise_default_plot_colors');
            plot(pp.xdatenums(1:nr)',d{:}(:,1),matlab_args{:})%,'color',[0,0,0]
            if nc>1
                % add some hair
                %---------------
                hold on
                iter=0;
                for tt=1:size(d{:},1)
                    iter=iter+1;
                    plot(pp.xdatenums(tt+(0:nc-1))',d{:}(tt,:)','color',map{iter})%,'linestyle','--'
                    if iter==numel(map)
                        iter=0;
                    end
                end
                hold off
            end
            set(gca,'xlim',pp.xlim,'XTick',pp.tickLocs,'XtickLabel',pp.xtick_labels,...
                'ylim',[ymin,max(ymax,ymin+eps)])
            vout={gca()};
        elseif strcmp(funct_type,{'plot_decomp'})
            uall=d{:};
            upos=zeros(size(uall));
            uneg=zeros(size(uall));
            ipos=uall>=0;
            ineg=uall<0;
            upos(ipos)=uall(ipos);
            uneg(ineg)=uall(ineg);
            b1=bar(pp.xdatenums,upos,'stack'); %  area(this_pos)
            hold on
            b2=bar(pp.xdatenums,uneg,'stack'); % area(this_neg)
            hold on
            for ii=1:numel(b1)
                set(b2(ii),'FaceColor',b1(ii).FaceColor)
            end
            plot(pp.xdatenums,sum(uall,2),'k-','linewidth',2);
            axis tight;
            hold off
%             set([b1,b2],'edgeColor','k')
            set(gca,'xlim',pp.xlim,'XTick',pp.tickLocs,'XtickLabel',pp.xtick_labels)
            vout={gca()};
        else % : plot, plotyy, bar, etc.
            if strcmp(funct_type,'plotyy')
                [ax12,h1,h2]=plotfunc(pp.xdatenums,d{1},pp.xdatenums,d{2},matlab_args{:});
                ax12={ax12};
                if ~isempty(linewidth)
                    set([h1,h2(:).'],linewidth{:})
                end
                 h1h2={h1,h2(:).'};
            else
               plotfunc(pp.xdatenums,d{:},matlab_args{:});
                ax12={gca};
            end
            hold on
            dall=[d{:}];
            dmin=utils.stat.nanmin(dall(:));
            dmax=utils.stat.nanmax(dall(:));
            add_vertical_lines();
            add_horizontal_lines();
            hold off
            try
                set(cell2mat(ax12),'xlim',pp.xlim,'XTick',pp.tickLocs,'XtickLabel',pp.xtick_labels)
            catch
                % 2014B vagaries
                set([ax12{:}],'xlim',pp.xlim,'XTick',pp.tickLocs,'XtickLabel',pp.xtick_labels)
            end
            if strcmp(funct_type,'plotyy')
                vout=[ax12,h1h2(1:nout-1)];
            else
                h=get(ax12{1},'Children');
                vout={h};
            end
        end
        varargout=vout(1:nout);
        
        function add_horizontal_lines()
            if ~isempty(rise_items.hline)
                xmin=pp.xlim(1);
                xmax=pp.xlim(2);
                for iline=1:numel(rise_items.hline)
                    hline=rise_items.hline(iline);
                    if ~(dmin<=hline && dmax>=hline)
                        warning(['horizontal line "',num2str(hline),'" outside the range of data plotted'])
                    end
                    plot([xmin,xmax],[hline,hline],'r',linewidth{:})
                end
            end
        end
        
        function add_vertical_lines()
            if ~isempty(rise_items.vline)
                for iline=1:numel(rise_items.vline)
                    vline=rise_items.vline(iline);
                    if pp.xlim(1)<=vline && pp.xlim(2)>=vline
                        plot([vline,vline],[dmin,dmax],'r',linewidth{:})
                    else
                        warning(['date "',char(serial2date(vline)),'" outside the range of data plotted'])
                    end
                end
            end
        end
    end

end

function [datta,date_numbers,sizdata]=stretch_data(datta,date_numbers,sizdata)

n=numel(datta);

if n==1
    
    return
    
end

lowestdn=inf;

highestdn=-inf;

for its=1:n
    
    lowestdn=min(lowestdn,date_numbers{its}(1));
    
    highestdn=max(highestdn,date_numbers{its}(end));
    
end

dn=lowestdn:highestdn;

T=numel(dn);

for its=1:n
    
    start=find(date_numbers{its}(1)==dn,1,'first');
    
    final=find(date_numbers{its}(end)==dn,1,'last');
    
    % override the first dimension
    sizdata{its}(1)=T;
    
    newdata=nan(sizdata{its});
    
    newdata(start:final,:,:,:)=datta{its};
    
    datta{its}=newdata;
    
    date_numbers{its}=dn;
    
end
    
end
