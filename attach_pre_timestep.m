%  Add user-specific calculations to be performed before each timestep.
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
%  Add user-specific calculations to be performed after each timestep.
%  This is not a function but rather lines of code that are inserted in 
%  the main routine at the right place. All variables known in the THERMAID.m
%  function are also known in this script.

%  If no additional pre-timestep calculations are necessary, this file can be empty.
%  However, it must exist in order for Thermaid to function properly.