function varargout=myplot(plotfunc,varargin)

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

pp=plot_specs(date_numbers{1},nticks,rise_items.date_format);

if rise_items.subplots
    [varargout{1:nargout}]=do_multiple_plots();
else
   [varargout{1:nargout}]=plot_it(datta);
end
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
            plotfunc(d{:},this.varnames,matlab_args{:});
            vout={gca()};
        elseif strcmp(funct_type,{'hist'})
            plotfunc(d{:},matlab_args{:});
            vout={gca()};
        else
            h1h2={};
            if strcmp(funct_type,{'plotyy'})
                [ax12,h1,h2]=plotfunc(pp.xdatenums,d{1},pp.xdatenums,d{2},matlab_args{:});
                h1h2={h1,h2};
                ax12={ax12};
                if ~isempty(linewidth)
                    set([h1,h2],linewidth{:})
                end
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
            set(cell2mat(ax12),'xlim',pp.xlim,'XTick',pp.tickLocs,'XtickLabel',pp.xtick_labels)
            vout=[ax12,h1h2(1:nout-1)];
        end
        %,'XTickMode','auto'...
        %     'ylim',[ymin,ymax])
        grid on
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

    function graph_handle=do_multiple_plots()
        r0=rise_items.figsize(1);
        c0=rise_items.figsize(2);
        nstar=r0*c0;
        npar=size(datta{1},2);
        nfig=ceil(npar/nstar);
        graph_handle=nan(nfig,1);
        dd=cell(1,ndatasets);
        for fig=1:nfig
            if nfig>1
                titelfig=[rise_items.figtitle,' ',int2str(fig)];
            else
                titelfig=rise_items.figtitle;
            end
            graph_handle(fig)=figure('name',titelfig);
            [Remains,r,c]=number_of_rows_and_columns_in_figure(fig,npar,r0,c0);
            for idata=1:min(nstar,Remains)
                subplot(r,c,idata)
                var_id=(fig-1)*nstar+idata;
                for ids=1:ndatasets
                    dd{ids}=datta{ids}(:,var_id);
                end
                plot_it(dd);
                if ndatasets==1
                    title(this{1}.varnames{var_id})
                end
            end
        end
    end
end
