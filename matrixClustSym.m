function p = matrixClustSym(varargin)
% MATRIXCLUSTSYM - Iterative clustering of a symmetric adjacency matrix
% matrixClustSym clusters a symmetric adjacency matrix by iteratively
% computing the eigenvector centrality --interpreting it as a cluster
% weight-- and using it to truncate the adjacency matrix for the next
% step. 
%
% USAGE:
%  p = matrixClustSym(s)
%  p = matrixClustSym(s, del)
%  p = matrixClustSym(s, del)
%  p = matrixClustSym(s, del, esz)
%  p = matrixClustSym(s, del, esz, out)
% 
%    s     - NxN sparse adjacency matrix
%    del   - cutoff, stop iterating if max(max(s))<=del (default 1e-2)
%    esz   - if network component has fewer than esz nodes, use eig
%            instead of eigs (default 1000)
%    out   - display info while running (out=1) or not (out=0) (default 1)
%    p     - matrix of cluster weights
%
% Written by
%   Tom Michoel
%   tom.michoel@frias.uni-freiburg.de

s = varargin{1};
switch nargin
 case 1
  del = 1e-2;
  esz = 1000;
  out = 1;
 case 2
  del = varargin{2};
  esz = 1000;
  out = 1;
 case 3
  del = varargin{2};
  esz = varargin{3};
  out = 1;
 case 4
  del = varargin{2};
  esz = varargin{3};
  out = varargin{4};
end

% some params
opteigs.disp = 0;
dim = length(s);

% we need nonzeros on the diagonal to get the disconnected components, see
% <http://blogs.mathworks.com/steve/2007/03/20/connected-component-labeling-part-3/>
s2 = s - diag(diag(s)); % first remove any possible elements on diagonal, then put ones everywhere
if issparse(s) 
    s2 = s2+speye(length(s2));
else
    s2 = s2+eye(length(s2));
end

% intialize probabilities and 1-weights to all-zeros
p = sparse(dim,1);
m = p;

% loop until all nodes are assigned with total probability at least 1-del
while max(s2(:))>del
  % create truncated adjacency matrix
  [I,J,sv] = find(s2);
  tw = sqrt(1-m);
  sv = tw(I).*tw(J).*sv;
  s2 = sparse(I,J,sv,dim,dim);
  % decompose s2 in disconnected components
  [a,b,c,d] = dmperm(s2);
  if length(c)==2 % one component
    if dim > esz 
      [v,ev,flag] = eigs(s2,1,'lm',opteigs);
      if flag~=0
        disp('eigs did not converge');
      end
    else
      [v,ev] = eig(full(s2));
      v = v(:,end);
      ev = ev(end,end);
    end
  else % multiple components
    s2per = s2(a,b); % block diagonally permuted s2
    ev = 0;
    ixmax = 0;
    for ix = 1:length(c)-1
      % block matrix
      blk = s2per(c(ix):c(ix+1)-1,d(ix):d(ix+1)-1);
      % 1st and last blocks don't correspond to separate component if not
      % square 
      if size(blk,1)==size(blk,2)
        if length(blk) > esz 
          [vtmp,etmp,flag] = eigs(blk,1,'lm',opteigs);
          if flag~=0
            disp('eigs did not converge');
          end
        else
          [vtmp,etmp] = eig(full(blk));
          vtmp = vtmp(:,end);
          etmp = etmp(end,end);
        end
        % check if we increase emax
        if etmp > ev
          ev = etmp;
          % embed block vector in whole space
          v = zeros(dim,1);
          v(c(ix):c(ix+1)-1) = vtmp;
          % permute indices back to original order
          [~,ainv] = sort(a);
          v = v(ainv);
        end
      end
    end
    if out
      disp(['CLUSTER ', num2str(size(p,2)), ', max(s)=', num2str(max(max(s2))), ...
            ', max(ev)=', num2str(ev), ', min(p)=', num2str(min(sum(p,2)))]);
    end
  end
  % create v with positive entries
  [~,imax] = max(abs(v));
  % normalize
  v = v./v(imax); 
  % this one has all positive elements, upto numerical precision, hence:
  v(v<0) = 0;
  % next column of p, largest weight should be scaled according to
  % total probability already assigned to imax
  p = [p min((1-m(imax)).*v,1-m)];
  % update weights
  m=m+p(:,end);
end
% first column of p is empty
p = p(:,2:end);
