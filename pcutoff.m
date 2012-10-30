function [q,n1,n2,nc] = pcutoff(p,x)
% PCUTOFF - set elements <=x to 0 in sparse matrix p
% [q,n1,n2,nc] = pcutoff(p,x)
%   p  - large sparse matrix
%   x  - cutoff
%   q  - sparse matrix returned
%   n1 - number of rows with at least one non-zero
%   n2 - number of rows with at least two non-zeros
%   nc - number of columns with more than one non-zero

  [a,b,pv] = find(p);
  c = find(pv>x);
  q = sparse(a(c),b(c),pv(c),size(p,1),size(p,2));
  q2 = spones(q);
  n1 = length(find(sum(q2,2)));
  n2 = length(find(sum(q2,2)>1));
  nc = length(find(sum(q2,1)>1));