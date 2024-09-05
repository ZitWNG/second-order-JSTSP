function [new_L, new_G] = attackNodeG(i, Ad, alpha, T)
[N,~] = size(Ad);
new_Ad = Ad;
new_Ad(i,:) = 0;    % 收不到消息
% new_Ad(:,i) = 0;    % 无法发送消息
indegree = sum(new_Ad, 2);  % 求入度（A的每行行和）
new_L = diag(indegree) - new_Ad;         % Laplacian 矩阵 L = D - A
G_A = eye(N) - T^2/2 * new_L;
G_B = T * eye(N) - alpha * T^2/2 * new_L;
G_C = - T * new_L;
G_D = eye(N) - alpha * T * new_L;
new_G = [G_A, G_B; G_C, G_D];
end