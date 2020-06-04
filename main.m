clean all;
close all;
clc;


model.n = 100;

model.T = 200;
n_sim = 40;

% Quality definition
model.q      = rand(model.n,1);
model.q      = sort(model.q, 'descend');

final_indegree=zeros(model.n, n_sim);
final_outdegree=zeros(model.n, n_sim);

for sim=1:n_sim
    A = simulate(model);
    final_indegree(:, sim) = sum(A,1)';
    final_outdegree(:, sim) = sum(A,2)';
end

%% Indegree Analysis

indegree_analysis(model, final_indegree);



%% Outdegree Analysis

outdegree_analysis(model, final_outdegree);


% fashion_blogger=[15.8E6, 7.6E6, 5.1E6, 5.1E6, 5.2E6, 2.7E6, 2E6, 1.7E6, 1.3E6, 1E6, 1E6]';
% n=fashion_blogger(1,1);
% m= size(fashion_blogger,1);
% expected_indegree= (n+1)*ones(m,1) - [n:-1:n-m+1]';
% expected_indegree= n./expected_indegree;
% 
% figure();
% scatter(fashion_blogger, [1:m]);
% hold on;
% plot(expected_indegree, [1:m], 'g', 'LineWidth',3);
% set(gca, 'XScale', 'log');
% set(gca, 'YScale', 'log');



% 
% 
% Z = ones(n, n) - X;
% Y = zeros(n, n+1);
% Y(:,1) = ones(n,1);
% Y(:,n+1) = ones(n,1);
% for d=0:n-1
%     for j=1:n
%         Y(d+1,1) = Y(d+1,1)*Z(j,d+1);
%         Y(d+1,n+1) = Y(d+1,n+1)*X(j,d+1);
%     end
% end
% for k=1:floor(n/2)
%     combos = nchoosek(1:n,k);
%     for com=1:size(combos,1)
%         temp_1 = ones(n,1);
%         temp_2 = ones(n,1);
%         for el=1:n
%             if (any(combos(com,:) == el))
%                 temp_1 = temp_1.*X(el,:)';
%                 temp_2 = temp_2.*Z(el,:)';
%             else
%                 temp_1 = temp_1.*Z(el,:)';
%                 temp_2 = temp_2.*X(el,:)';
%             end
%         end
%         Y(:,k+1) = Y(:,k+1) + temp_1;
%         if (n-k ~= k)
%             Y(:,n-k+1) = Y(:,n-k+1) + temp_2;
%         end
%     end
% end
% 
% expected_numberOfNodes = Y*linspace(0,n,n+1)';
% expected_frequency    = expected_numberOfNodes/n;
% 
% Y = Y/n;
% 
% [D,K]=meshgrid(0:1:n-1, 0:1:n);
% figure();
% surf(D,K,Y');
% shading interp;
% hold on;
% view(2);
% ylim([0,n]);
% xlim([1,n-1]);
% set(gca, 'YScale', 'log');
% set(gca, 'XScale', 'log');
% ylabel('frequency');
% xlabel('indegree');
% hold off;

% figure();
% bins = xbins(2:end-1);
% %xbins = [xbins(1,1), xbins(1,3), (xbins(3:end-1)+xbins(4:end))/2];
% xbins = [xbins(1,1), (xbins(3:end-1)+xbins(4:end))/2];
% a=histcounts(indegree(:,T), xbins);
% %a=histcounts(indegree(:,T));
% % for ii=1:size(b,1)
% %     c(ii)=sum(a(ii:end));
% % end
% % a = c;
% 
% %scatter(b, a/sum(a), 'filled');
% %hold on;
% % plot(linspace(0,n-1,n)', expected_frequency, 'r');
% f = fit(bins(3:end)', a(2:end)'/sum(a), 'power1')
% plot(f,bins(3:end)', a(2:end)'/sum(a));
% % f = fit(linspace(1,n-2,n-2)', expected_frequency(2:n-1), 'power1');
% % plot(f,linspace(1,n-2,n-2)', expected_frequency(2:n-1));
% set(gca, 'YScale', 'log');
% set(gca, 'XScale', 'log');
% ylabel('frequency');
% xlabel('indegree');
% hold off;
% 
% figure();
% hist(indegree(:,T));
% 
% 
% 
% figure();
% hist(outdegree(:,T));
% hold on;
% outdeg = outdegree_prob(n)*n;
% plot(outdeg(:,1), 'r');
% hold off;
% 
% toc


