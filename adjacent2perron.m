function P = adjacent2perron(A)
[N,~] = size(A);
indegree = sum(A,2);        % 求入度（A的每行行和）
L = diag(indegree) - A;     % Laplacian 矩阵 L = D - A
ep = 1/(5*max(indegree));	% 该值越小 收敛越慢数据越多
P = eye(N,N) - ep * L;      % Perron 矩阵 P = I - eL
end