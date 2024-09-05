clear
clc;
close all
warning('off');

%% 图的生成和初始状态设计
%%% 3个结点固定
N = 3;
Ad = [0 1 0; 1 0 1; 0 1 0];
init_px = [600 600 600]'; init_py = [1000 1600 2400]';
targ_px = [600 1200 600]'; targ_py = [1200 1600 2000]';
init_vx = [300 300 0]'; init_vy = [0 -300 300]';

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
init_x = [init_px - targ_px; init_vx];
init_y = [init_py - targ_py; init_vy];
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
phi = 1;%0.99;

%% 普通迭代
global_x1 = init_x;
global_y1 = init_y;
[n,~] = size(G);
N = n / 2;
for k = 1:k_limits
    global_x1(:,k+1) = G * global_x1(:,k);
    global_y1(:,k+1) = G * global_y1(:,k);
end
[output, L_final] = estimation_2(global_x1, L, T, alpha);
disp(L_final);
disp("x1");
disp(output);
% [xxx, L_final] = estimation_2(global_x1, L, T, 2);
% disp(xxx);

%% 机制加噪1
global_x2 = init_x;
global_y2 = init_y;
[n,~] = size(G);
N = n/2;
b = rand(N, k_limits);      % start indicator 
theta_x2 = zeros(N, k_limits);
theta_y2 = zeros(N, k_limits);
Tc = 2;
amp = -100;
for k = 1:k_limits
    for i = 1:N
        if k < k_limits - 10 && b(i, k) > 0.5
            theta_x2(i, k) = amp * phi^k + theta_x2(i, k);
            theta_y2(i, k) = amp * phi^k + theta_y2(i, k);
            % Tc = randi([2 5]);          % compensation period
            theta_x2(i, k + Tc) = amp * phi^k + theta_x2(i, k + Tc);
            theta_y2(i, k + Tc) = amp * phi^k + theta_y2(i, k + Tc);
            for l = 1:Tc - 1
                theta_x2(i, k + l) = -2 * amp * phi^k / (Tc - 1) + theta_x2(i, k + l);
                theta_y2(i, k + l) = -2 * amp * phi^k / (Tc - 1) + theta_y2(i, k + l);
            end
        end
    end
    input_x2 = [T^2/2 * theta_x2(:, k); T * theta_x2(:, k)];
    input_y2 = [T^2/2 * theta_y2(:, k); T * theta_y2(:, k)];
    global_x2(:,k+1) = G * global_x2(:,k) + input_x2;
    global_y2(:,k+1) = G * global_y2(:,k) + input_y2;
end
[output2, L_final2] = estimation_2(global_x2, L, T, alpha);
% disp(L_final2);
disp("x2");
disp(output2);

%% 机制加噪2
global_x3 = init_x;
global_y3 = init_y;
[n,~] = size(G);
N = n/2;
b = rand(N, k_limits);      % start indicator
theta_x3 = zeros(N, k_limits);
theta_y3 = zeros(N, k_limits);
amp2 = 100;
for k = 1:k_limits
    for i = 1:N
        if k < k_limits - 10 && b(i, k) > 0.1 %k < k_limits / 2 && b(i, k) > 0.9
            theta_x3(i, k) = amp2 * phi^k + theta_x3(i, k);
            theta_y3(i, k) = amp2 * phi^k + theta_y3(i, k);
            Tm = 1;%randi([1 5]);  % compensation period
            Tn = 3;%randi([Tm + 1 Tm + 5]);
            theta_x3(i, k + Tm) = -Tn / (Tn-Tm) * amp2 * phi^k + theta_x3(i, k + Tm);
            theta_x3(i, k + Tn) = Tm / (Tn-Tm) * amp2 * phi^k + theta_x3(i, k + Tn);
            theta_y3(i, k + Tm) = -Tn / (Tn-Tm) * amp2 * phi^k + theta_y3(i, k + Tm);
            theta_y3(i, k + Tn) = Tm / (Tn-Tm) * amp2 * phi^k + theta_y3(i, k + Tn);
        end
    end
    input_x3 = [T^2/2 * theta_x3(:, k); T * theta_x3(:, k)];
    input_y3 = [T^2/2 * theta_y3(:, k); T * theta_y3(:, k)];
    global_x3(:,k+1) = G * global_x3(:,k) + input_x3;
    global_y3(:,k+1) = G * global_y3(:,k) + input_y3;
end
[output3, L_final3] = estimation_2(global_x3, L, T, alpha);
% disp(L_final3);
disp("x3");
disp(output3);

%% 绘图
% %%% 静态 tvx 图
% figure;
% for i=1:N
%     plot(1:k_limits,global_x1(i,1:end-1),'LineWidth',1.5);xlim([0 k_limits]);ylim([-400 400]);hold on;grid on;
%     set(gca,'FontSize',16);
%     xlabel('Time-step','FontSize',16);ylabel('Positions of the Agents','FontSize',16);
% end
%%% 静态 tvy 图
figure;
for i = 1:N
    plot(1:k_limits,global_x1(N + i,1:end-1),'LineWidth',2);xlim([0 k_limits]);ylim([0 400]);hold on;grid on;
    set(gca,'FontSize',16);
    xlabel('Time-step','FontSize',16);ylabel('Velocities of the Agents','FontSize',16);
end
%%% 静态 xy 图
figure;
for i = 1:N
    plot(global_x1(i,1:end-1) + targ_px(i) * ones(1, k_limits), global_y1(i,1:end-1) + targ_py(i) * ones(1, k_limits), 'LineWidth',2);xlim([0 3500]);ylim([500 2500]);hold on;grid on;
    set(gca,'FontSize',16);
    xlabel('X Axis','FontSize',16);ylabel('Y Axis','FontSize',16);
end
plot([global_x1(1,1) + targ_px(1), global_x1(2,1) + targ_px(2)], [global_y1(1,1) + targ_py(1), global_y1(2,1) + targ_py(2)], '-o', 'Color', 'b', 'LineWidth', 2);
plot([global_x1(3,1) + targ_px(3), global_x1(2,1) + targ_px(2)], [global_y1(3,1) + targ_py(3), global_y1(2,1) + targ_py(2)], '-o', 'Color', 'b', 'LineWidth', 2);
plot([global_x1(1,1) + targ_px(1), global_x1(3,1) + targ_px(3)], [global_y1(1,1) + targ_py(1), global_y1(3,1) + targ_py(3)], '-o', 'Color', 'b', 'LineWidth', 2);
plot([global_x1(1,end-1) + targ_px(1), global_x1(2,end-1) + targ_px(2)], [global_y1(1,end-1) + targ_py(1), global_y1(2,end-1) + targ_py(2)], '-o', 'Color', 'b', 'LineWidth', 2);
plot([global_x1(3,end-1) + targ_px(3), global_x1(2,end-1) + targ_px(2)], [global_y1(3,end-1) + targ_py(3), global_y1(2,end-1) + targ_py(2)], '-o', 'Color', 'b', 'LineWidth', 2);
plot([global_x1(1,end-1) + targ_px(1), global_x1(3,end-1) + targ_px(3)], [global_y1(1,end-1) + targ_py(1), global_y1(3,end-1) + targ_py(3)], '-o', 'Color', 'b', 'LineWidth', 2);
% %%% 动态 xy 图
% figure;
% color = [[0.00,0.45,0.74];[0.85,0.33,0.10];[0.93,0.69,0.13];[0.72,0.27,1];[0.47,0.67,0.19]];
% for k = 1:k_limits
%     for i = 1:N
%         plot([global_x1(i,k)+ targ_px(i),global_x1(i,k+1)+ targ_px(i)],[global_y1(i,k)+ targ_py(i),global_y1(i,k+1)+ targ_py(i)],'-.','MarkerSize',10,'Color',color(i,:));xlim([0 5000]);ylim([0 3000]);
%         hold on;
%     end
% %     pause(0.05);   % 暂停，就可以看到点的变化走向
% end
%%% 静态 tvy 图
figure;
for i = 1:N
    plot(1:k_limits,global_x2(N + i,1:end-1),'LineWidth',2);xlim([0 k_limits]);ylim([0 400]);hold on;grid on;
    set(gca,'FontSize',16);
    xlabel('Time-step','FontSize',16);ylabel('Velocities of the Agents','FontSize',16);
end

%%% 静态 xy 图
figure;
for i = 1:N
    plot(global_x2(i,1:end-1) + targ_px(i) * ones(1, k_limits), global_y2(i,1:end-1) + targ_py(i) * ones(1, k_limits), 'LineWidth',2);xlim([0 3500]);ylim([500 2500]);hold on;grid on;
    set(gca,'FontSize',16);
    xlabel('X Axis','FontSize',16);ylabel('Y Axis','FontSize',16);
end
plot([global_x2(1,1) + targ_px(1), global_x2(2,1) + targ_px(2)], [global_y2(1,1) + targ_py(1), global_y2(2,1) + targ_py(2)], '-o', 'Color', 'b', 'LineWidth', 2);
plot([global_x2(3,1) + targ_px(3), global_x2(2,1) + targ_px(2)], [global_y2(3,1) + targ_py(3), global_y2(2,1) + targ_py(2)], '-o', 'Color', 'b', 'LineWidth', 2);
plot([global_x2(1,1) + targ_px(1), global_x2(3,1) + targ_px(3)], [global_y2(1,1) + targ_py(1), global_y2(3,1) + targ_py(3)], '-o', 'Color', 'b', 'LineWidth', 2);
plot([global_x2(1,end-1) + targ_px(1), global_x2(2,end-1) + targ_px(2)], [global_y2(1,end-1) + targ_py(1), global_y2(2,end-1) + targ_py(2)], '-o', 'Color', 'b', 'LineWidth', 2);
plot([global_x2(3,end-1) + targ_px(3), global_x2(2,end-1) + targ_px(2)], [global_y2(3,end-1) + targ_py(3), global_y2(2,end-1) + targ_py(2)], '-o', 'Color', 'b', 'LineWidth', 2);
plot([global_x2(1,end-1) + targ_px(1), global_x2(3,end-1) + targ_px(3)], [global_y2(1,end-1) + targ_py(1), global_y2(3,end-1) + targ_py(3)], '-o', 'Color', 'b', 'LineWidth', 2);

% %%% 静态 xy 图
% figure;
% for i = 1:N
%     plot(global_x3(i,1:end-1) + targ_px(i) * ones(1, k_limits), global_y3(i,1:end-1) + targ_py(i) * ones(1, k_limits), 'LineWidth',1.5);xlim([0 5000]);ylim([0 3000]);hold on;grid on;
%     set(gca,'FontSize',16);
%     xlabel('Time-step','FontSize',16);ylabel('Velocities of the Agents','FontSize',16);
% end


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



