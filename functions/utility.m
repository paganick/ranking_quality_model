function utility = utility(A, i, q, hub, c, d)
    %utility = A(i,:)*q;
%    utility = 1000*max(A(i,:)'.*q);



    utility = max(A(i,:)'.*q);
    
%     followees = A(i,:)'.*q;
%     followees(followees == 0) = [];
%     if (size(followees,1) > 0)
%         utility = sum(followees)/size(followees,1);
%     else
%         utility = 0;
%     end



%     utility =  utility + d*diversity(A,i);
%     if (i == hub)
%         q = update_q(A, q, hub);
%         %utility = utility - c*cost(A,i); % + d*diversity(A,i);
%         utility = q(i) - c*cost(A,i);
%     else
%         utility = utility - c*(sum(A(i,:)))^2; % + d*diversity(A,i);
%     end    
end

