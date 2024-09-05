clear;
close all

%%
N = 3;
GN = 2 * N;
Ad = [0 1 1; 0 0 4; 0 2 0]; % [0 1 0; 0 0 1; 1 0 0]; % [0 0 0; 1 0 1; 0 0 0];% 需要是一个有 spanning tree 的有向图
init_px = [800 1200 800]'; init_py = [1000 1600 2400]';
targ_px = [600 1200 600]'; targ_py = [1200 1600 2000]';
init_vx = [300 300 300]'; init_vy = [300 100 150]';

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
k_limits = 10000;


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
phi = 1;%0.99;

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
    Y = global_y2(:,1:k);
    conY(1,k) = cond(Y*Y');
    yn6(1,k) = norm(inv(Y*Y'),2);
    yn7(1,k) = norm((Y),2);%norm((Y*Y'),2);
    yn8(1,k) = norm(Y'*inv(Y*Y'),2);
    tmp = sort(eig(Y*Y'));
    y85(1,k) = tmp(1);
    yn9(1,k) = norm(input_x2(:,1:k),2);
    ynn(1,k) = norm(input_x2(:,1:k) * Y'*inv(Y*Y'),2);
end
figure;
plot(20:k_limits, conY(1,20:end),'r');hold on
plot(20:k_limits, 12*(20:k_limits).^(2.3),'b');

% figure;
% plot(20:k_limits, yn6(1,20:end),'r');hold on
% plot(20:k_limits, 12*(20:k_limits).^(-1),'b');

% figure;
% plot(20:k_limits, yn7(1,20:end),'r');hold on
% plot(20:k_limits, 12*(20:k_limits).^(2.95),'b');

figure;
plot(20:k_limits, yn8(1,20:end),'r');hold on
plot(20:k_limits, 8*(20:k_limits).^(-0.55),'b');

figure;
plot(20:k_limits, y85(1,20:end),'r');hold on
% plot(20:k_limits, 8*(20:k_limits).^(0.5),'b');

figure;
plot(20:k_limits, yn9(1,20:end),'r');hold on
plot(20:k_limits, 8*(20:k_limits).^(0.5),'b');

figure;
plot(20:k_limits, ynn(1,20:end),'r');hold on
plot(20:k_limits, 18*(20:k_limits).^(0.05),'b');

