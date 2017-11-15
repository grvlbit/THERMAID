function str = time2str(time)
%  Converts time in seconds to a human readable string
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

if time<86400 % Simulation < 1 day
    str = [num2str(time) ' s'];
elseif time < 365*86400 % Simulation < 1 year
    days = floor(time/86400);
    seconds = time - 86400*days;
    str = [num2str(days) ' d ' num2str(seconds) ' s' ];
else
    years = floor(time/(365*86400));
    tmp = time - years*365*86400;
    days = floor(tmp/86400);
    seconds = tmp - 86400*days;
    str = [num2str(years) ' a ' num2str(days) ' d ' num2str(seconds) ' s'];
end
end