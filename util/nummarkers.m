%------------------------BEGIN FUNCTION
% NUMMARKERS takes a vector of line handles in h
% and reduces the number of plot markers on the lines
% to num. This is useful for closely sampled data.
%
% example:
% t = 0:0.01:pi;
% p = plot(t,sin(t),'-*',t,cos(t),'r-o');
% nummarkers(p,10);
% legend(p,'sin(t)','cos(t)')
%
% Magnus Sundberg Feb 08, 2001
%
% Note : if you plot in a for loop, store the p=plot(..) into an array plVec, and pass this to this function, calling it once only!
%h = arrayof plot handles
%num = numer of desired markers
%doSpacingAlongX = 1 (default) -> marker delta-x constant, 0-> spacing constant along the curve
%
% Massimo Ciacci, Nov 04, 2011

function nummarkers(h,num, doSpacingAlongX)
if nargin < 3
    doSpacingAlongX = 1;
end

for n = 1:length(h)
    if strcmp(get(h(n),'type'),'line')
        axes(get(h(n),'parent'));
        x = get(h(n),'xdata');
        y = get(h(n),'ydata');
        if (doSpacingAlongX)
            ti = round(linspace(1,length(x),num));
        else
            t = 1:length(x);
            s = [0 cumsum(sqrt(diff(x).^2+diff(y).^2))]; %measures length along the curve
            si = (0:num-1)*s(end)/(num-1); %equally spaced lengths along the curve
            si(end) = s(end); %fix last point to be within the curve
            ti = round(interp1(s,t,si)); %find x index of markers
        end
        xi = x(ti);
        yi = y(ti);
        marker = get(h(n),'marker');
        markerS = get(h(n),'MarkerSize');
        color = get(h(n),'color');
        style = get(h(n),'linestyle');
        width = get(h(n),'linewidth');
        % make a line with just the markers
        set(line(xi,yi),'marker',marker,'MarkerSize',markerS,'linestyle','none','color',color,'linewidth',width);
        % make a copy of the old line with no markers
        set(line(x,y),'marker','none','linestyle',style,'linewidth',width,'color',color);
        % set the x- and ydata of the old line to [], this tricks legend to keep on working
        set(h(n),'xdata',[],'ydata',[]);
    end
end