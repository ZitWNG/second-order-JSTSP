function Ar = perron2adjacent(P)
% [N,~] = size(P);
% epL = eye(N,N) - P;
% Ar = zeros(N);
% for x = 1 : N-1 % 不知道\epsilon 具体是多少，所以试一下最大入度是多少
%     L = epL * 5 * x;
%     if max(max(L)) - 1 > x && max(max(L)) + 1 < x
%         continue
%     else
%         [row,col] = find(-1.5 < L & L < -0.5);
%         [m,~] = size(row);
%         for i = 1:m
%             Ar(row(i),col(i)) = 1;
%         end
%     end
% end

%% new
% 其实还有很多别的可以考虑的，入度和Ar的1匹配
[N,~] = size(P);
epL = eye(N,N) - P;
Ar = zeros(N);
preL = epL * 5;
minvu = N;
for x = 1:N
    if preL(x,x) < minvu
        minvu = preL(x,x);
    end
end
L = preL.*round(1/minvu);

for i = 1:N
    for j = 1:N
        pairsum = L(i,j) + L(j,i);
        if pairsum < -2/N % && pairsum >-1.5
            Ar(i,j) = 1;
            Ar(j,i) = 1;
        end
    end
end

end