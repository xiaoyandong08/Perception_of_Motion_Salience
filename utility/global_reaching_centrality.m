function [Gr,Cr] = global_reaching_centrality(A)

d = distances(digraph(A));
d(1:size(d,1)+1:size(d,1)^2) = inf;
Cr = sum(1./d,2)./(size(d,1)-1);
Gr = sum(max(Cr) - Cr)/(size(d,1)-1);

end