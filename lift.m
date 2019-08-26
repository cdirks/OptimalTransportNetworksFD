%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Functional lifting                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function one_u = lift(u,range_u)
[m,n] = size(u);
one_u = zeros(m,n,length(range_u));
for i = 1:m
    for j = 1:n
        for k = 1:length(range_u)
            if (u(i,j) > range_u(k) )             
                one_u(i,j,k) = 1;
            end
        end
    end
end
end