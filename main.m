clear
clc;
close all
warning('off');

%% 图的生成和初始状态设计
%%% 3个结点固定
N = 3;
Ad = [0 1 0; 1 0 1; 0 1 0];
init_p = [2 4 6]';
targ_p = [6 3 2]';
init_v = [1 0 1]';

%%% 5个节点固定
% N = 5;
% Ad = [0 2 0 0 2; 2 0 3 1 3; 0 3 0 1 0; 0 1 1 0 2; 2 3 0 2 0];
% init_p = [-26 -12 15 31 -8]';   % 初始状态 0
% targ_p = [-26 -3 13 28 17]';   % 初始状态 5
% init_v = [1 0 2 3 -2]';

%%% 随机生成无向图
% N = 10;
% Ad = ERrandomGraph(N,0.6);
% init_p = round((rand(N,1)-0.5)*100);   % 初始状态 0
% targ_p = round((rand(N,1)-0.5)*100);   % 初始状态 5
% init_v = round((rand(N,1)-0.5)*10);

%%% 随机生成有向图，内部可定义节点的度范围
% N = 10;
% A = time_varying_digraph(N);
% A = (A > 0 | A < 0); % 将A中不等于0的元素记为1
% init_state = round((rand(N,1)-0.5)*100); %% [-50 50];

%% Laplacian 矩阵构建
%%% 求入度（A的每行行和）
indegree = sum(Ad, 2);
L = diag(indegree) - Ad;         % Laplacian 矩阵 L = D - A

% x_c = calc_con_point(P, init_state);    % 求 consensus point

%% 迭代参数设置
% k_limits = N^2;
k_limits = 100;
repetition = 10;

%% 系统设置
T = 0.1; 
alpha = 2; 
init_x = [init_p - targ_p; init_v];
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
phi = 0.99;

%% 普通迭代
global_x1 = inputDesign_0(init_x, k_limits, G, T);
[output, L_final] = estimation_2(global_x1, L, T, alpha);
disp(L_final);
disp("x1");
disp(output);
% [xxx, L_final] = estimation_2(global_x1, L, T, 2);
% disp(xxx);

%% 机制加噪1
global_x2 = inputDesign_1(init_x, k_limits, G, T, phi, 2);
[output2, L_final2] = estimation_2(global_x2, L, T, alpha);
% disp(L_final2);
disp("x2");
disp(output2);

%% 机制加噪2
global_x3 = inputDesign_1(init_x, k_limits, G, T, phi, 3);
[output3, L_final3] = estimation_2(global_x3, L, T, alpha);
% disp(L_final3);
disp("x3");
disp(output3);

%% 绘图
figure;
for i=1:N
    plot(1:k_limits,global_x1(i,1:end-1),'LineWidth',2);xlim([0 k_limits]);ylim([-25 10]);hold on;grid on;
    set(gca,'FontSize',16);
    xlabel('Time-step','FontSize',16);ylabel('Positions of the Agents','FontSize',16);
end
figure;
for i = 1:N
    plot(1:k_limits,global_x1(N + i,1:end-1),'LineWidth',1.5);xlim([0 k_limits]);ylim([-10 20]);hold on;grid on;
    set(gca,'FontSize',16);
    xlabel('Time-step','FontSize',16);ylabel('Velocities of the Agents','FontSize',16);
end
figure;
for i = 1:N
    plot(1:k_limits,global_x2(i,1:end-1),'LineWidth',1.5);xlim([0 k_limits]);ylim([-25 10]);hold on;grid on;
    set(gca,'FontSize',16);
    xlabel('Time-step','FontSize',16);ylabel('Positions of the Agents','FontSize',16);
end
figure;
for i = 1:N
    plot(1:k_limits,global_x2(N + i,1:end-1),'LineWidth',1.5);xlim([0 k_limits]);ylim([-10 20]);hold on;grid on;
    set(gca,'FontSize',16);
    xlabel('Time-step','FontSize',16);ylabel('Velocities of the Agents','FontSize',16);
end
figure;
for i = 1:N
    plot(1:k_limits,global_x3(i,1:end-1),'LineWidth',1.5);xlim([0 k_limits]);ylim([-25 10]);hold on;grid on;
    set(gca,'FontSize',16);
    xlabel('Time-step','FontSize',16);ylabel('Positions of the Agents','FontSize',16);
end
figure;
for i = 1:N
    plot(1:k_limits,global_x3(N + i,1:end-1),'LineWidth',1.5);xlim([0 k_limits]);ylim([-10 20]);hold on;grid on;
    set(gca,'FontSize',16);
    xlabel('Time-step','FontSize',16);ylabel('Velocities of the Agents','FontSize',16);
end

% figure;
% for i=1:N
%     subplot(3,2,1);plot(1:k_limits,global_x1(i,1:end-1),'LineWidth',1.5);xlim([0 k_limits]);ylim([-25 10]);hold on;title('position original design');grid on;
%     subplot(3,2,2);plot(1:k_limits,global_x1(N + i,1:end-1),'LineWidth',1.5);xlim([0 k_limits]);ylim([-10 20]);hold on;title('velocity original design');grid on;
%     subplot(3,2,3);plot(1:k_limits,global_x2(i,1:end-1),'LineWidth',1.5);xlim([0 k_limits]);ylim([-25 10]);hold on;title('position input design 1');grid on;
%     subplot(3,2,4);plot(1:k_limits,global_x2(N + i,1:end-1),'LineWidth',1.5);xlim([0 k_limits]);ylim([-10 20]);hold on;title('velocity input design 1');grid on;
%     subplot(3,2,5);plot(1:k_limits,global_x3(i,1:end-1),'LineWidth',1.5);xlim([0 k_limits]);ylim([-25 10]);hold on;title('position input design 2');grid on;
%     subplot(3,2,6);plot(1:k_limits,global_x3(N + i,1:end-1),'LineWidth',1.5);xlim([0 k_limits]);ylim([-10 20]);hold on;title('velocity input design 2');grid on;
% end


% %% 画图
% %%% 画配色自定义的图
% b = bar(data,'FaceColor','flat');
% newc = [0.03,0.41,0.67;0.94,0.98,0.91;0.73,0.89,0.74;0.48,0.8,0.77;0.26,0.63,0.79];
% for k = 1:size(data,2)
%     b(k).CData = newc(mod(k,5)+1,:);
% end
% ylabel('Frobenius norm error performance','FontSize',18);
% set(gca, 'XTickLabel', {'\alpha = 0','\alpha = 0.3','\alpha = 0.6','\alpha = 0.9','\alpha = 1.2'} );
% set(gca,'FontSize',18)
% hold on
% 
% %%% draw the lines
% figure;
% for nn = 1:alpha_num
%     plot(dist_6(nn,:),output_6(nn,:,2),'LineWidth',1.5);hold on;
% end
% grid on;
% legend();



