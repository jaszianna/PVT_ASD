%%This code belogs to Balázs Hangya's lab
%%We used this code to plot the firing activity of the given neuron

function plot_raster_lines_fast(actTime,yCoords,lineColor)
%FAST_RASTERPLOT(ACTTIME) plots the activity of a cell, firing at the
%timepoints specified in ACTTIME vector. Improved performance because
%raster lines are not handled as individual objects but as parts of a
%continous line containing NaNs between each line.
%
%   Optional parameters:
%   YCOORDS: 2 element row vector, y coordinates of the raster lines (default: [0,1]).
%   LINECOLOR: colorcode (character) for plot color (default: k -> black).

if nargin<3
    if nargin<2
        yCoords = [0,1];
    end
    lineColor = 'k';
end

actTime = reshape(actTime,1,[]); %convert to row vector

xVector = [actTime;actTime;NaN(1,length(actTime))];
xVector = xVector(:);
yVector = repmat([yCoords.';NaN],1,length(actTime));
yVector = yVector(:);
plot(xVector,yVector,'Color',lineColor)
end