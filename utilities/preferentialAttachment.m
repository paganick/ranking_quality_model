function el = preferentialAttachment(n,m)
%PREFERENTIALATTACHMENT simple preferential attachment for network growth.

% The probability that a new vertex attaches to a given old vertex
%                           is proportional to the (total) vertex degree.

% @input n final number of desired vertices
% @input m # of edges to attach at each step
% @output el, edge list in Mx3 format. 

% Note 1: Vertices arrive one at a time.
% Note 2: Assume undirected simple graph.
% Source: Newman, "The Structure and Function of Complex Networks"
%         B-A., "Emergence of Scaling in Random Networks"

% Other routines used: weightedRandomSample.m
% Updated: fprintf, documentation

% IB: last updated, 3/24/14
%##################################################################



vertices = 2;
if not(vertices<=n); fprintf('Specify more than 2 nodes.\n');  return; end
el=[1 2 1; 2 1 1];      % start with one edge


while vertices < n
  
  vertices=vertices+1;  % add new vertex

  if m>=vertices
    for node=1:vertices-1
      el = [el; node vertices 1];
      el = [el; vertices node 1];  % add symmetric edge
    end
    continue
  end
    
  deg=[];        % compute nodal degrees for this iteration
  for v=1:vertices; deg=[deg; (numel(find(el(:,2)==v)))]; end
  

  % add m edges
  r = weightedRandomSample(m,[1:vertices],deg/sum(deg));
  while not(length(unique(r))==length(r))
    r = weightedRandomSample(m,[1:vertices],deg/sum(deg));
  end
  
  for node=1:length(r)
    el = [el; r(node) vertices 1];
    el = [el; vertices r(node) 1];      
  end      
  
end


function s = weightedRandomSample(n,P,W)
%WEIGHTEDRANDOMSAMPLE Weighted random sampling.
% 
% INPUTs: number of draws from a discrete distribution (n)
%         possible values to pick from, (P)
%         set of normalized weights/probabilities, (W)
% OUTPUTs: s - set of n numbers drawn from P
%              according to the weights in W


    s = [];

    if abs(sum(W)-1)>10^(-8); fprintf('The probabilities do not sum up to 1.\n'); return; end

    % divide the unit interval into |P| segments each with length W_i
    unit = [0,W(1)];
    for w=2:length(W)
      unit = [unit W(w)+unit(length(unit))]; %#ok<AGROW>
    end

    % draw a random number in the unit interval uniformly - where does it fall?
    while length(s)<n

      lb = find(unit<rand, 1, 'last' );    % rand is in [unit(lb), unit(lb+1)]
      s = [s P(lb)];             %#ok<AGROW> % pick P(lb)

    end
