function [x1,y1,x2,y2] = read_dfn_data_from_file(filename,format,c1,c2,c3,c4, startRow, endRow)
%  Read DFN data from file. Currently FracMan and FracSim3D CSV files are
%  supported. 
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
%  Authors: Gunnar Jansen, University of Neuchatel, 2016-2018
%
%  [x1,y1,x2,y2] = read_dfn_data_from_file(filename,format,c1,c2,c3,c4, startRow, endRow)
%
%  Input: 
%        filename          name of the file to read (string)
%        format            format of the input file (string that is either
%                                 'fracman' for FracMan files 
%                              or 'fracsim' for FracSim3D files
%        c1,c2,c3,c4       fracture position columns to be read from file
%                                  in order:
%                                      xStart
%                                      yStart
%                                      xEnd
%                                      yEnd
%        startRow          (optional) first row to be read from file
%        endRow            (optional) last row to be read from file
%
%  Output:
%        x1,y1,x2,y2       fracture position columns
%                                  in order:
%                                      xStart
%                                      yStart
%                                      xEnd
%                                      yEnd

%% Initialize variables.
delimiter = ',';
if nargin<=6
    startRow = 3;
    endRow = inf;
end

%% Read columns of data as strings:
% For more information, see the TEXTSCAN documentation.
if (strcmp(format,'fracman'))
    formatSpec = '%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%s';
elseif (strcmp(format,'fracsim'))
    formatSpec = '%f%f%f%f%f%f';
else
    error('Unknown input format')
end

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
% dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1);
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines', startRow(1)-1);

for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'HeaderLines', startRow(block)-1);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

x1 = cell2mat(dataArray(:, c1));
y1 = cell2mat(dataArray(:, c2));
x2 = cell2mat(dataArray(:, c3));
y2 = cell2mat(dataArray(:, c4));

%% Close the text file.
fclose(fileID);