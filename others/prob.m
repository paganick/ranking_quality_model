clear all;
close all;
clc;

N=200;
k=0;
tic
outdeg = outdegree_prob(N);
toc
figure();
plot(outdeg(:,1));
xlim([1, 20]);

% for k=1:N-1
%     for d=1:N-1
%         indeg(k, d) = indegree_prob(N,d,k);
%     end
% end
% [X,Y]=meshgrid(1:1:N-1);
% figure();
% surf(X,Y,indeg);


function p = outdegree_prob(N)
p=zeros(N-1,N-1);
p(1,:)=1/(N-1)*ones(1,N-1);
for d=2:N-1
    for k_0=1:N-1
        for k=k_0:N-d
            p(d,k_0) = p(d,k_0)+1/(N-1-k)*p(d-1,k+1);
        end
    end
end
end

function p = indegree_prob(N,d,k)
b = nchoosek(N-1, d);
p = b *1/(N-k)^d * (1-1/(N-k))^(N-1-d);
end


