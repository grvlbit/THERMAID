%  Add user-defined calculations or visualization commands to be performed
%  after each timestep.
%  ---------------------------------------------------------------------
%  Copyright (C) 2017 by the Thermaid authors
% 
%  This file is part of Thermaid.
% 
%  Thermaid is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
% 
%  Thermaid is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License
%  along with Thermaid.  If not, see <http://www.gnu.org/licenses/>.
%  ---------------------------------------------------------------------
% 
%  Authors: Gunnar Jansen, University of Neuchatel, 2017 
%
%  Add user-defined calculations or visualization commands to be performed
%  after each timestep. This is not a function but rather lines of code
%  that are inserted in the main routine at the right place. All variables
%  known in the THERMAID.m function are also known in this script.
%  This file includes the default runtime visialization.  

%  Generate figures for the pressure and temperature solutions with handles.
if showPlot && time == dt0
    s = linspace(0,1,256);
    rgb1=[0.230, 0.299, 0.754];
    rgb2=[0.706, 0.016, 0.150];
    cmap = diverging_map(s,rgb1,rgb2);

    figure(1)
    h1= pcolor(x,y,udata.T0');
    hold on
    colormap(cmap)
    line([XY1(:,1)';XY1(:,3)'],[XY1(:,2)';XY1(:,4)'],'Color','r');
    h6 = scatter(xe,ye,40,pf,'filled','s', 'MarkerEdgeColor',[.1 .1 .1],...
                'LineWidth',0.5);
    shading interp
    axis equal
    xlabel('x [m]')
    ylabel('y [m]')
    title('Pressure')
    
    if(udata.flagHeatTransport)
        figure(3)
        h4 = pcolor(x,y,tNew');
        hold on
        colormap(cmap)
        shading interp
        axis equal
        caxis([0 udata.tmax]); 
        line([XY1(:,1)';XY1(:,3)'],[XY1(:,2)';XY1(:,4)'],'Color','k','LineWidth',2);
        h5 = scatter(xe,ye,30,tNewf,'filled');
        drawnow
    end
end

% Update the data shown in the figures based on their handles. The axis and labels remain unchanged. This is much faster than plotting the whole figure each time.
if showPlot
        set(h1,'cdata',p');
        set(h6,'cdata',pf);
        
        if(udata.flagHeatTransport)
            set(h4,'cdata',tNew');
            set(h5,'cdata',tNewf);
        end
        drawnow
end

% Optional data output to file after each timestep
%
%     outfile = ['data_' num2str(time)];
%     save(['run/isc/output/' outfile],'p','pf','cNew','cNewf','tNew','tNewf','trig_tog')

