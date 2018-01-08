function [row,col] = intersectionsSegments(XY1,XY2)
%  Finds the intersections of the fracture segments with each other and
%  return the row and column indices of the interesctions
%  ---------------------------------------------------------------------
%  Copyright (C) 2016 by the Thermaid authors
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
%  Authors: Gunnar Jansen, University of Neuchatel, 2016-2017
%
%  intersectionsSegments(XY1,XY2)
%
%  Input:
%        XY1      (nf, 4)       1st array containing all fracture segments as
%                               [x_begin y_begin x_end y_end]
%        XY2      (nf, 4)       2nd array containing all fracture segments as
%                               [x_begin y_begin x_end y_end]
%
%  Output:
%        row    (#ffIntersections,1) row indices of fracture-fracture
%                               intersections
%        col    (#ffIntersections,1) col indices of fracture-fracture
%                               intersections
%
%  Usage:
%        In Thermaid the two input arrays are identical to find the
%        intersections of the fracture set with itself.


row = [];
col = [];

for i=1:length(XY1)
    for j=1:length(XY2)
        s1_x = XY1(i,3) - XY1(i,1);     s1_y = XY1(i,4) - XY1(i,2);
        s2_x = XY1(j,3) - XY1(j,1);     s2_y = XY1(j,4) - XY1(j,2);
        
        if ((-s2_x * s1_y + s1_x * s2_y) ~= 0)
            s = (-s1_y * (XY1(i,1) - XY1(j,1)) + s1_x * (XY1(i,2) - XY1(j,2))) / (-s2_x * s1_y + s1_x * s2_y);
            t = ( s2_x * (XY1(i,2) - XY1(j,2)) - s2_y * (XY1(i,1) - XY1(j,1))) / (-s2_x * s1_y + s1_x * s2_y);
        else
            s = Inf;
            t = Inf;
        end
        
        if (s > 0 && s < 1 && t > 0 && t < 1) %if (s >= 0 && s <= 1 && t >= 0 && t <= 1) % The commented version included merely touching segments.
            touch = true;
        else
            touch = false; % No collision
        end
        if (touch)
            row = [row i j];
            col = [col j i];
        end
    end
end

end

