clear;
close all
% phi = 1;
% nn = zeros(1,1000);
% for r = 1:1:100
%     x = ones(3,1);
%     for k = 1:1:1000
%         if mod(k,3) == 1
%             add = phi^k;
%         elseif mod(k,3) == 2
%             add = -2 * phi^k;
%         else
%             add = phi^k;
%         end
%         x(:, end+1) = x(:, end) + ones(3, 1) + add;
%         nn(1,k) = nn(1,k) + norm(x'*pinv(x*x'),'fro');
%     end
% end
% figure;
% plot(nn,'b');hold on
% x = 1:1:1000;
% y = 1./sqrt(x);
% plot(x, y,'r');

%%
% for k = 1:100
%     A1 = 1:2:2*k;
%     A2 = ones(1,k);
%     A = [A1, A2];
%     yn1(1,k) = norm(A*A','fro');
%     yn2(1,k) = norm(inv(A*A'),'fro');
%     yn3(1,k) = yn1(1,k) * yn2(1,k) * k^2.5;
%     yn4(1,k) = cond(A*A');      % condition number
%     yn5(1,k) = norm(A*A',2);
%     yn6(1,k) = norm(inv(A*A'),2);
%     yn7(1,k) = norm(A*A',2) * norm(inv(A*A'),2);
% end
% figure;
% plot(1:100,yn1,'b');hold on
% plot(1:100,yn2,'r');hold on
% plot(1:100,yn3,'g');
% % figure;
% % plot(1:100,yn5,'b');hold on
% % plot(1:100,yn6,'r');hold on
% % plot(1:100,yn7,'g');

%%
N = 3;
GN = 2 * N;
Ad = [0 1 1; 0 0 4; 0 2 0]; % [0 1 0; 0 0 1; 1 0 0]; % [0 0 0; 1 0 1; 0 0 0];% 需要是一个有 spanning tree 的有向图
init_px = [600 1200 600]'; init_py = [1200 1600 2000]';
targ_px = [600 1200 600]'; targ_py = [1200 1600 2000]';
init_vx = [300 300 300]'; init_vy = [300 300 300]';

%%% Laplacian 矩阵构建
indegree = sum(Ad, 2);  % 求入度（A的每行行和）
L = diag(indegree) - Ad;         % Laplacian 矩阵 L = D - A
%%% alpha、T 设置和验证
T = 0.1;
alpha = 2;
%%% G 矩阵构建
G_A = eye(N) - T^2/2 * L;
G_B = T * eye(N) - alpha * T^2/2 * L;
G_C = - T * L;
G_D = eye(N) - alpha * T * L;
G = [G_A, G_B; G_C, G_D];
%%% 系统初始化
init_x = [init_px - targ_px; init_vx];
init_y = [init_py - targ_py; init_vy];

%%% 迭代参数设置
k_limits = 2000;

%%
% global_x1 = init_x;
% global_y1 = init_y;
% 
% for k = 1: k_limits
%     global_x1(:,k+1) = G * global_x1(:,k);
%     global_y1(:,k+1) = G * global_y1(:,k);
% end
% [output, L_final] = estimation_2(global_x1, L, T, alpha);
% [outputG, G_final] = estimation(global_x1, G);
% L_calcu = (eye(N) - G_final(1:N, 1:N)) * 2 / (T^2);
% for k = 1:k_limits
%     Y = T * global_x1(1:N,1:k) + alpha * T * global_x1(N+1:2*N,1:k) + rand(N, k) - 0.5 * ones(N, k);
%     conY(1,k) = cond(Y*Y');
%     yn6(1,k) = norm(inv(Y*Y'),2);
%     yn7(1,k) = norm((Y*Y'),2);
%     yn8(1,k) = norm(Y'*inv(Y*Y'),2);
% end
% figure;
% % plot(20:k_limits, conY(1,20:end),'r');hold on
% % plot(20:k_limits, 12*(20:k_limits).^(2.3),'b');
% 
% % plot(20:k_limits, yn6(1,20:end),'r');hold on
% % plot(20:k_limits, 12*(20:k_limits).^(-1),'b');
% 
% % plot(20:k_limits, yn7(1,20:end),'r');hold on
% % plot(20:k_limits, 12*(20:k_limits).^(2.95),'b');
% 
% plot(20:k_limits, yn8(1,20:end),'r');hold on
% plot(20:k_limits, 12*(20:k_limits).^(-0.8),'b');
% 
% % plot(20:k_limits, yn7(1,20:end).*yn6(1,20:end)./conY(1,20:end),'r');hold on
% % plot(20:k_limits, 0.1*(20:k_limits).^(0),'b');
% % % plot(20:k_limits, yn6(1,20:end),'b');


%%
global_x2 = init_x;
global_y2 = init_y;
b = rand(N, k_limits);      % start indicator 
theta_x2 = zeros(N, k_limits);
theta_y2 = zeros(N, k_limits);
a_x2 = zeros(N, k_limits);
a_y2 = zeros(N, k_limits);
Tc = 2;
amp = 100;
phi = 0.99;%1;

for k = 1:k_limits
    for i = 1:N
        if k < k_limits - 40 && b(i, k) > 0.5
            a_x2(i, k) = amp * phi^k * rand(1);
            a_y2(i, k) = amp * phi^k * rand(1);
            theta_x2(i, k) = a_x2(i, k) + theta_x2(i, k);
            theta_y2(i, k) = a_y2(i, k) + theta_y2(i, k);
            % Tc = randi([2 5]);          % compensation period
            for l = 1:Tc - 1
                theta_x2(i, k + l) = -2 * a_x2(i, k) / (Tc - 1) + theta_x2(i, k + l);
                theta_y2(i, k + l) = -2 * a_y2(i, k) / (Tc - 1) + theta_y2(i, k + l);
            end
            theta_x2(i, k + Tc) = a_x2(i, k) + theta_x2(i, k + Tc);
            theta_y2(i, k + Tc) = a_y2(i, k) + theta_y2(i, k + Tc);
        end
    end
    attacked = 1;
    input_x2(:, k) = [T^2/2 * theta_x2(:, k); T * theta_x2(:, k)];
    input_y2(:, k) = [T^2/2 * theta_y2(:, k); T * theta_y2(:, k)];

    global_x2(:,k+1) = G * global_x2(:,k) + input_x2(:, k);
    global_y2(:,k+1) = G * global_y2(:,k) + input_y2(:, k);
end
for k = 1:k_limits
    Y = T * global_x2(1:N,1:k) + alpha * T * global_x2(N+1:2*N,1:k);
    conY(1,k) = cond(Y*Y');
    yn6(1,k) = norm(inv(Y*Y'),2);
    yn7(1,k) = norm((Y*Y'),2);
    yn8(1,k) = norm(Y'*inv(Y*Y'),2);
    tmp = sort(eig(Y*Y'));
    y85(1,k) = tmp(1);          % 最小特征值
    [U,S,V] = svd(Y);
    y86(1,k) = norm(input_x2(N+1:2*N,1:k)*(V)');  % \theta V'
    y87(1,k) = 1/sqrt(y85(1,k)) * y86(1,k);
    yn9(1,k) = norm(input_x2(N+1:2*N,1:k),2);
    ynn(1,k) = norm(input_x2(N+1:2*N,1:k) * Y'*inv(Y*Y'),2);
end
% figure;
% plot(20:k_limits, conY(1,20:end),'r');hold on
% plot(20:k_limits, 12*(20:k_limits).^(2.3),'b');
% figure;
% plot(20:k_limits, yn6(1,20:end),'r');hold on
% plot(20:k_limits, 12*(20:k_limits).^(-1),'b');
% figure;
% plot(20:k_limits, yn7(1,20:end),'r');hold on
% plot(20:k_limits, 12*(20:k_limits).^(2.95),'b');
% figure;
% plot(20:k_limits, yn8(1,20:end),'r');hold on
% plot(20:k_limits, 2*(20:k_limits).^(-0.55),'b');

figure;
plot(20:k_limits, y85(1,20:end),'r');hold on
% plot(20:k_limits, 8*(20:k_limits).^(0.5),'b');

figure;
plot(20:k_limits, y86(1,20:end),'r');hold on
plot(20:k_limits, 8*(20:k_limits).^(0.5),'b');

figure;
plot(20:k_limits, y87(1,20:end),'r');hold on
plot(20:k_limits, 8*(20:k_limits).^(0.5),'b');

% figure;
% plot(20:k_limits, yn9(1,20:end),'r');hold on
% % plot(20:k_limits, 8*(20:k_limits).^(0.5),'b');

figure;
plot(20:k_limits, ynn(1,20:end),'r');hold on
plot(20:k_limits, 18*(20:k_limits).^(0.05),'b');

