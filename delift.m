%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Functional delifting                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u = delift(v,range)
[k,l,~] = size(v);
hz = range(2)-range(1);
% Number of pixels equal to 1 in 3rd dimension
u_pixels = zeros(k,l);
u = zeros(k,l);
for i = 1:k
    for j = 1:l
        u_pixels(i,j) = sum(v(i,j,:));
        u(i,j) = range(1) + u_pixels(i,j)*hz;
    end
end
end
