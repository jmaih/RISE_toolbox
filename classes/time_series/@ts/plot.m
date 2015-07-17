function plot_handle=plot(varargin)
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
% More About
% ------------
%
% Examples
% ---------
%
% See also: 


%     plot(X,db) plots time series db versus vector of dates X. X can be of
%     multiple forms and need not be of the same length as db:
%         - a string of dates. e.g. X='1990:2000' X='1990Q1:2000Q1'
%         X='1990m1:2000m2' 
%         - serial dates. e.g. X=date2serial('1990Q1):date2serial('2000Q1')
%         X='1990m1:2000m2' 
%         - If X is a scalar, disconnected line objects are created and
%           plotted as discrete points vertically at X.
%  
%     plot(db) plots the columns of db versus the dates within db
%     If db is complex, plot(db) is equivalent to plot(real(db),imag(db)).
%     In all other uses of plot, the imaginary part is ignored.
%  
%     Various line types, plot symbols and colors may be obtained with
%     plot(X,db,S) where S is a character string made from one element
%     from any or all the following 3 columns:
%  
%            b     blue          .     point              -     solid
%            g     green         o     circle             :     dotted
%            r     red           x     x-mark             -.    dashdot 
%            c     cyan          +     plus               --    dashed   
%            m     magenta       *     star             (none)  no line
%            y     yellow        s     square
%            k     black         d     diamond
%            w     white         v     triangle (down)
%                                ^     triangle (up)
%                                <     triangle (left)
%                                >     triangle (right)
%                                p     pentagram
%                                h     hexagram
%                           
%     For example, plot(X,db,'c+:') plots a cyan dotted line with a plus 
%     at each data point; plot(X,db,'bd') plots blue diamond at each data 
%     point but does not draw any line.
%  
%     plot(X1,Y1,S1,X2,Y2,S2,X3,Y3,S3,...) combines the plots defined by
%     the (X,db,S) triples, where the X's and db's are vectors or matrices 
%     and the S's are strings.  
%  
%     For example, plot(X,db,'y-',X,db,'go') plots the data twice, with a
%     solid yellow line interpolating green circles at the data points.
%  
%     The plot command, if no color is specified, makes automatic use of
%     the colors specified by the axes ColorOrder property.  By default,
%     plot cycles through the colors in the ColorOrder property.  For
%     monochrome systems, plot cycles over the axes LineStyleOrder property.
%  
%     Note that RGB colors in the ColorOrder property may differ from
%     similarly-named colors in the (X,db,S) triples.  For example, the 
%     second axes ColorOrder property is medium green with RGB [0 .5 0],
%     while plot(X,db,'g') plots a green line with RGB [0 1 0].
%  
%     If you do not specify a marker type, plot uses no marker. 
%     If you do not specify a line style, plot uses a solid line.
%  
%     plot returns a column vector of handles to lineseries objects, one
%     handle per plotted line. 
%  
%     The X,db pairs, or X,db,S triples, can be followed by 
%     parameter/value pairs to specify additional properties 
%     of the lines. For example, plot(X,db,'LineWidth',2,'Color',[.6 0 0]) 
%     will create a plot with a dark red line width of 2 points.
%     
%     In addition to those matlab properties, RISE adds further properties,
%     which allow to control for. See parse_plot_args

% plot('1994m7:1997m1',db(:,:,1),...
%     'figsize',[2,2],...
%     'figtitle','no title',...
%     'nticks',10,...
%     'legend',{'v1','v2'},...
%     'legend_loc','BO',...
%     'logy',true,...
%     'secondary_y',{'v1','v4'},...
%     'subplots',true,...
%     'linewidth',2,...
%     'date_format',17);
% xrotate(90)

plot_handle0=utils.plot.myplot(@plot,varargin{:});
if nargout
    plot_handle=plot_handle0;
end

end