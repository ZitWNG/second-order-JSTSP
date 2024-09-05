clear
clc;
close all
warning('off');

%% 图的生成和初始状态设计
%%% 3个结点固定
Ad = [0 1 0; 1 0 1; 0 1 0];
init_p = [2 4 6]';
init_v = [1 0 0]';
targ_p = [6 3 2]';

%% Laplacian 矩阵构建
%%% 求入度（A的每行行和）
indegree = sum(Ad, 2);
L = diag(indegree) - Ad;         % Laplacian 矩阵 L = D - A

%% 迭代参数设置
% k_limits = N^2;
k_limits = 200;
repetition = 100;

%% 系统设置
T = 0.1;
alpha = 2;
init_x = [init_p - targ_p; init_v];
[N, ~] = size(Ad);
G_A = eye(N) - T^2/2 * L;
G_B = T * eye(N) - alpha * T^2/2 * L;
G_C = - T * L;
G_D = eye(N) - alpha * T * L;
G = [G_A, G_B; G_C, G_D];
A = [1, T; 0, 1];
BK = [T^2/2, alpha * T^2/2; T, alpha * T];
lambda = eig(L);
% 计算 asymptotic convergence factor
rasym = 0;
for i = 1:N-1
    rasym = max(rasym, max(abs(eig(A-lambda(i+1) * BK))));
end

%% 普通迭代
global_x1 = init_x;
for k = 1:k_limits
    global_x1(:,k+1) = G * global_x1(:,k);
end

[output, L_final] = estimation_2(global_x1, L, T, alpha);
disp(output);

%% 机制加噪1

ansp = zeros(1, 50);

for r = 1:repetition
    j = 1;
    global_x2 = init_x;
    b = rand(N, k_limits);      % start indicator
    % phi = 0.99;                 % attenuation factor
    theta = zeros(N, k_limits);
    for phi = 0.5:0.01:0.99
        for k = 1:k_limits
            for i = 1:N
                if mod(k, 10) == 1
                    theta(i, k) = phi^k + theta(i, k);
                    Tc = randi([2 5]);          % compensation period
                    theta(i, k + Tc) = phi^k + theta(i, k + Tc);
                    for l = 1:Tc - 1
                        theta(i, k + l) = -2 * phi^k / (Tc - 1) + theta(i, k + l);
                    end
                end
            end
            input = [T^2/2 * theta(:, k); T * theta(:, k)];
            global_x2(:,k+1) = G * global_x2(:,k) + input;
        end
        [output2, Les] = estimation_2(global_x2, L, T, alpha);
        ansp(1,j) = ansp(1,j) + output2(2)/repetition;
        j = j + 1;
    end
end
plot(0.5:0.01:0.99,ansp(1,1:end),'-.');
% disp(output2);
% i = 1;
% anss = zeros(1, 21);
% outt = zeros(1, 21);
% for x = -1:0.1:1
%     theta = [-x;-x;x];
%     input = [T^2/2 * theta; T * theta];
%     global_x2(:,k_limits +2) = G * global_x2(:,k_limits +1) + input;
%     global_x2(:,k_limits +3) = G * global_x2(:,k_limits +2) - 2 * input;
%     global_x2(:,k_limits +4) = G * global_x2(:,k_limits +3) + input;
%     [output2, ~] = estimation_2(global_x2, L, T, alpha);
%     anss(1,i) = output2(2);
%     outt(1,i) = calc_dist(global_x2);
%     i = i + 1;
%
% end
% figure;
% plot(-1:0.1:1,anss(1,1:end),'-.');
% figure;
% plot(1:21,outt(1,1:end));

%% 求逆
x_1a = global_x2(1:N, 1:end-1);
x_1b = global_x2(1:N, 2:end);
x_2a = global_x2(N+1:2*N, 1:end-1);
x_2b = global_x2(N+1:2*N, 2:end);
Za_1 = T^2/2 * x_1a + alpha * T^2/2 * x_2a;
Za_2 = T * x_1a + alpha * T * x_2a;
Zb_1 = x_1a + T * x_2a - x_1b;
Zb_2 = x_2a - x_2b;
Za = [Za_1, Za_2];
Zb = [Zb_1, Zb_2];
L = Zb * Za' * inv(Za * Za');
% disp(L);

%% 约束求逆
AA = ones(N, 1);
xx = sum(sum(inv(Za * Za')));
xy = [1;1;1] * xx * [1 1 1];
L2 = Zb * Za' * inv(Za * Za') * (eye(N) - inv(Za * Za') * AA * inv(AA' * inv(Za * Za') * AA) * AA')';
% disp(L2);

%% 渐进
global_x4 = init_x;
sigma = 1;
fi = 0.99;
theta4 = zeros(N, k_limits);
diss = zeros(1, repetition);
outt = zeros(1, repetition);
for r = 1:repetition
    rand_list = sigma * randn(N, k_limits);  % 高斯分布
    for k = 1:k_limits
        if k == 1
            theta4(:,k) = rand_list(:, k);
        else
            theta4(:,k) = fi^k * rand_list(:, k) - fi^(k-1) * rand_list(:, k-1);
        end
        input4 = [T^2/2 * theta4(:, k); T * theta4(:, k)];
        global_x4(:,k+1) = G * global_x4(:,k) + input4;
    end
    [output4, L4] = estimation_2(global_x4, L, T, alpha);
    diss(1,r) = calc_dist(global_x4);
    outt(1,r) = output4(2);
end
% figure;
% scatter(outt, diss);

%% 机制加噪2
% global_x3 = init_x;
% b = rand(N, k_limits);      % start indicator
% phi = 0.99;                 % attenuation factor
% theta = zeros(N, k_limits);
% for k = 1:k_limits
%     for i = 1:N
%         if k < k_limits / 2 && b(i, k) > 0.9
%             theta(i, k) = phi^k + theta(i, k);
%             Tm = randi([1 5]);  % compensation period
%             Tn = randi([Tm + 1 Tm + 5]);
%             theta(i, k + Tm) = -Tn / (Tn-Tm) * phi^k + theta(i, k + Tm);
%             theta(i, k + Tn) = Tm / (Tn-Tm) * phi^k + theta(i, k + Tn);
%         end
%     end
%     input = [T^2/2 * theta(:, k); T * theta(:, k)];
%     global_x3(:,k+1) = G * global_x3(:,k) + input;
% end
%
% [output3, L_final3] = estimation_2(global_x3, L, T, alpha);
% disp(output3);
%
% %% 绘图
% figure;
% for i=1:N
%     subplot(3,2,1);plot(1:k_limits,global_x1(i,1:end-1));hold on;title('position x1');grid on;
%     subplot(3,2,2);plot(1:k_limits,global_x1(N + i,1:end-1));hold on;title('velocity x1');grid on;
%     subplot(3,2,3);plot(1:k_limits,global_x2(i,1:end-1));hold on;title('position x2');grid on;
%     subplot(3,2,4);plot(1:k_limits,global_x2(N + i,1:end-1));hold on;title('velocity x2');grid on;
%     subplot(3,2,5);plot(1:k_limits,global_x3(i,1:end-1));hold on;title('position x3');grid on;
%     subplot(3,2,6);plot(1:k_limits,global_x3(N + i,1:end-1));hold on;title('velocity x3');grid on;
% end





