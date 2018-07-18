function [G,d,channel_indexes] = makeGForConcavity(c)

% Read channel points:

[from_ind, from_val, to_ind, to_val, area] = textread('G_channels.txt', '%f, %f, %f, %f, %f');

% Read outlet points:

[outlet_index, elevation] = textread('G_outlets.txt','%f, %f');

% Allocate sparse matrix:

max_index = max([from_ind;to_ind;outlet_index])+1;
G = sparse(max_index,max_index);
d = zeros(length(outlet_index),1);

% Build sparse matrix for channel points:

for i=1:length(from_ind)
   G(i,from_ind(i)+1) = from_val(i) .* area(i).^c;
   G(i,to_ind(i)+1) = to_val(i) .* area(i).^c;
   channel_indexes(i) = from_ind(i)+1;
end

% Build sparse matrix for outlets:

for i=1:length(outlet_index)
    G(i+length(from_ind),outlet_index(i)+1) = 1.0;
    d(i) = elevation(i);
    outlet_indexes(i) = outlet_index(i)+1;
end

