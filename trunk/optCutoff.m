function [xopt, qopt, n1, n2, nc] = optCutoff(p,s)
% OPTCUTOFF - find optimal probability cutoff in fuzzy clusters
% [xopt, qopt] = optCutoff(p,s)
%   p    - sparse matrix of fuzzy clusters
%   s    - sparse adjacency matrix clustered by p
%   xopt - optimal cutoff for each cluster
%   qopt - p with elements below xopt set to 0, if xopt=0 whole column of
%          p is set to 0
%   n1   - number of rows with at least one non-zero in qopt
%   n2   - number of rows with at least two non-zeros in qopt
%   nc   - number of columns with at least one non-zero in qopt
 
xopt = zeros(1,size(p,2));
q = zeros(size(p));
for k=1:size(p,2)
  % cutoff values to test
  x = floor(1e3.*p(:,k))./1e3;
  x = unique(x(x>0));
  % find maximum
  [fmax,imax] = max(score(x,p(:,k),s));
  xopt(k) = x(imax);
  if xopt(k)>0
    q(:,k) = pcutoff(p(:,k),xopt(k));
  else
    q(:,k) = zeros(size(p(:,k)));
  end
end

% assign nodes to only one cluster
qopt = zeros(size(q));
[C,I] = max(q,[],2);
for k=1:size(q,1)
    qopt(k,I(k)) = C(k);
end
%qopt = q;

q2 = spones(qopt);
n1 = length(find(sum(q2,2)));
n2 = length(find(sum(q2,2)>1));
nc = length(find(sum(q2,1)>1));


function y = score(x,p,s)
% SCORE - score function to optimize cutoff 
%  score computes the ratio of the 'number of edges' at a given cutoff to
%  the 'number of nodes'; the score is optimal for the largest fully
%  connected cluster. 
  
[I,J,sv] = find(s);
y = zeros(size(x));
for k=1:length(x)
  ip = find(p>x(k));
  if ~isempty(ip)
    q = zeros(size(p));
    q(ip) = p(ip);
    y(k) = sum(q(I).*q(J).*sv)/sum(q);
  end
end