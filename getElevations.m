function [elevations,i] = getElevations()

[i, elevations] = textread('G_channel_elevations.txt','%f, %f');
