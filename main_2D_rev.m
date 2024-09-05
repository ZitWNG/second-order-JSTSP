clear
clc;
close all
warning('off');

%% 图的生成和初始状态设计
%%% 3个结点固定
N = 3;
GN = 2 * N;
Ad = [0 1 1; 0 0 4; 0 2 0]; % [0 1 0; 0 0 1; 1 0 0]; % [0 0 0; 1 0 1; 0 0 0];% 需要是一个有 spanning tree 的有向图
init_px = [600 600 600]'; init_py = [1000 1600 2400]';
targ_px = [600 1200 600]'; targ_py = [1200 1600 2000]';
init_vx = [300 300 0]'; init_vy = [-100 -200 100]';

%%% 5个节点固定
% N = 5;
% Ad = [0 2 0 0 2; 2 0 3 1 3; 0 3 0 1 0; 0 1 1 0 2; 2 3 0 2 0];
% init_p = [-26 -12 15 31 -8]';   % 初始状态 0
% targ_p = [-26 -3 13 28 17]';   % 初始状态 5
% init_v = [1 0 2 3 -2]';


%% 系统设置
%%% Laplacian 矩阵构建
indegree = sum(Ad, 2);  % 求入度（A的每行行和）
L = diag(indegree) - Ad;         % Laplacian 矩阵 L = D - A
%%% alpha、T 设置和验证
T = 0.1;
alpha = 2;
if test_alphaT(L, alpha, T) == 1
    disp("ok");
end
%%% G 矩阵构建
G_A = eye(N) - T^2/2 * L;
G_B = T * eye(N) - alpha * T^2/2 * L;
G_C = - T * L;
G_D = eye(N) - alpha * T * L;
G = [G_A, G_B; G_C, G_D];
%%% 系统初始化
init_x = [init_px - targ_px; init_vx];
init_y = [init_py - targ_py; init_vy];
%%% 计算 asymptotic convergence factor
eigenvalues = sort(eig(G));
rasym = eigenvalues(end-2);
%%% 迭代参数设置
k_limits = 200;    % N^2;
repetition = 10;
%%% 收敛状态计算
[con_vx, con_vy, con_px, con_py] = calcu_consensus(L, init_vx, init_vy, init_px, init_py, targ_px, targ_py, k_limits, T);

%% 普通迭代
k_attack =  k_limits; %20; %10;
[global_x1, global_y1] = inputDesign_ori(init_x, init_y, k_attack, G);
% [output, L_final] = estimation_2(global_x1, L, T, alpha);
[outputG, G_final] = estimation(global_x1, G);
L_calcu = (eye(N) - G_final(1:N, 1:N)) * 2 / (T^2);
print1 = ["x1:",outputG];
disp(print1);

%%% 找最重要结点
mmax = -inf;
for i = 1:N
    if L_calcu(i,i)> mmax
        mmax = L_calcu(i,i);
        node = i;   
    end
end
%%% 发起攻击
[new_L, new_G] = attackNodeG(node, Ad, alpha, T);   % 对结点发起攻击
% [new_L, new_G] = attackLinkG(2, 3, Ad, alpha, T);   % 对链路发起攻击
for k = k_attack + 1:k_limits
    global_x1(:,k+1) = new_G * global_x1(:,k);  %global_x1(node+N,k+1) = global_x1(node+N,k+1) + 1;% global_x1(node+N, k+1) = 0;
    global_y1(:,k+1) = new_G * global_y1(:,k);  %global_y1(node+N,k+1) = global_y1(node+N,k+1) + 1;% global_y1(node+N, k+1) = 0;
end
%%% distortion 和 deviation 分别表示 formation 的形变和 consensus state 的偏移
[distortion, deviation] = attack_performance(con_vx, con_vy, con_px, con_py, global_x1, global_y1, targ_px, targ_py, alpha);
printP = ["distortion:",distortion,"deviation:",deviation];
disp(printP);

% trav_w = [zeros(1,N),ones(1,1)] * pinv([L,ones(N,1)]);
% for samplek = 5:5:k_limits
%     con_px = (trav_w * (init_px - targ_px) + (samplek - 1) * T * trav_w * init_vx) * ones(N, 1) + targ_px;
%     con_py = (trav_w * (init_py - targ_py) + (samplek - 1) * T * trav_w * init_vy) * ones(N, 1) + targ_py;
%     last_px = global_x1(1:N,samplek);last_py = global_y1(1:N,samplek);                            % 最后时刻的 relative position
%     last_vx = global_x1(N+1:end,samplek);last_vy = global_y1(N+1:end,samplek);                    % 最后时刻的 velocity
%     delta_p = sqrt((con_px - targ_px - last_px).^2 + (con_py - targ_py - last_py).^2);  % 与 consensus position 的差异
%     delta_v = sqrt((con_vx - last_vx).^2 + (con_vy - last_vy).^2);                      % 与 consensus velocity 的差异
%     for i = 1:N
%         deviation(i, samplek/5) = delta_p(i,1) + alpha * delta_v(i,1);
%     end
% end
% figure;
% for i = 1:N
%     plot(5:5:k_limits, deviation(i, :), 'LineWidth',2);hold on;grid on;
% %     xlim([0 3500]);ylim([800 2500]);
%     set(gca,'FontSize',16);
%     xlabel('Time-step','FontSize',16);ylabel('Deviation','FontSize',16);
% end
% legend('Agent 1','Agent 2','Agent 3');

%% 机制加噪
seq2 = 1*[1 -1 -1 1];
[global_x2, global_y2] = inputDesign_noisy(init_x, init_y, k_attack, G, T, seq2);
[outputG2, G_final2] = estimation(global_x2, G);
L_calcu2 = (eye(N) - G_final2(1:N, 1:N)) * 2 / (T^2);

%%% 找最重要结点
mmax2 = -inf;
for i = 1:N
    if L_calcu2(i,i)> mmax2
        mmax2 = L_calcu2(i,i);
        node2 = i;   
    end
end
%%% 发起攻击
[new_L2, new_G2] = attackNodeG(node2, Ad, alpha, T);   % 对结点发起攻击
% [new_L2, new_G2] = attackLinkG(1, 3, Ad, alpha, T);   % 对链路发起攻击
for k = k_attack + 1:k_limits
    global_x2(:,k+1) = new_G2 * global_x2(:,k); % global_x2(node2+N,k+1) = global_x2(node2+N,k+1) + 1; % global_x2(node2+N, k+1) = 0;
    global_y2(:,k+1) = new_G2 * global_y2(:,k); % global_y2(node2+N,k+1) = global_y2(node2+N,k+1) + 1; % global_y2(node2+N, k+1) = 0; %
end
%%% distortion 和 deviation 分别表示 formation 的形变和 consensus state 的偏移
[distortion2, deviation2] = attack_performance(con_vx, con_vy, con_px, con_py, global_x2, global_y2, targ_px, targ_py, alpha);
printP2 = ["distortion:",distortion2,"deviation:",deviation2];
disp(printP2);

%%% deviation 图
% for samplek = 5:5:k_limits
%     con_px = (trav_w * (init_px - targ_px) + (samplek - 1) * T * trav_w * init_vx) * ones(N, 1) + targ_px;
%     con_py = (trav_w * (init_py - targ_py) + (samplek - 1) * T * trav_w * init_vy) * ones(N, 1) + targ_py;
%     last_px = global_x2(1:N,samplek);last_py = global_y2(1:N,samplek);                            % 最后时刻的 relative position
%     last_vx = global_x2(N+1:end,samplek);last_vy = global_y2(N+1:end,samplek);                    % 最后时刻的 velocity
%     delta_p = sqrt((con_px - targ_px - last_px).^2 + (con_py - targ_py - last_py).^2);  % 与 consensus position 的差异
%     delta_v = sqrt((con_vx - last_vx).^2 + (con_vy - last_vy).^2);                      % 与 consensus velocity 的差异
%     for i = 1:N
%         deviation2(i, samplek/5) = delta_p(i,1) + alpha * delta_v(i,1);
%     end
% end
% figure;
% for i = 1:N
%     plot(5:5:k_limits, deviation2(i, :), 'LineWidth',2);hold on;grid on;
% %     xlim([0 3500]);ylim([800 2500]);
%     set(gca,'FontSize',16);
%     xlabel('Time-step','FontSize',16);ylabel('Deviation','FontSize',16);
% end
% legend('Agent 1','Agent 2','Agent 3');

% seq3 = [1 -2 1];
% [global_x3, global_y3] = inputDesign_noisy(init_x, init_y, k_limits, G, T, seq3);
% seq4 = [1 0 -3 2];
% [global_x4, global_y4] = inputDesign_noisy(init_x, init_y, k_limits, G, T, seq4);
% [output2, L_final2] = estimation_2(global_x2, L, T, alpha);
% [outputG2, G_final2] = estimation(global_x2, G);

% L_calcu2 = (eye(N) - G_final2(1:N, 1:N)) * 2 / (T^2);
% ann = norm(G_final2-G,'fro')/norm(L_calcu2-L,'fro');
% disp(ann);
% print2 = ["x2:",output2];
% disp(print2);

% for k = 20:k_limits
%     [outputG2, G_final2] = estimation(global_x2(:,1:k), G);
%     L_calcu2 = (eye(N) - G_final2(1:N, 1:N)) * 2 / (T^2);
%     ans2(1,k) = norm(L_calcu2 - L,'fro');
%     
%     [outputG3, G_final3] = estimation(global_x3(:,1:k), G);
%     L_calcu3 = (eye(N) - G_final3(1:N, 1:N)) * 2 / (T^2);
%     ans3(1,k) = norm(L_calcu3 - L,'fro');
%     
%     [outputG4, G_final4] = estimation(global_x4(:,1:k), G);
%     L_calcu4 = (eye(N) - G_final4(1:N, 1:N)) * 2 / (T^2);
%     ans4(1,k) = norm(L_calcu4 - L,'fro');
% end
% figure;
% plot(20:k_limits, ans2(1,20:end),'b','LineWidth',2); hold on
% plot(20:k_limits, ans3(1,20:end),'r','LineWidth',2); hold on
% plot(20:k_limits, ans4(1,20:end),'g','LineWidth',2); hold on

% LLL_final2(1:N, 1:N) = (eye(N) - G_final2(1:N, 1:N))*2/T^2;
% LLL_final2(1:N, N+1:2*N) = (T * eye(N) - G_final2(1:N, N+1:2*N))*2/T^2;
% LLL_final2(N+1:2*N, N+1:2*N) = (eye(N) - G_final2(N+1:2*N, N+1:2*N))/T;
% LLL_final2(N+1:2*N, 1:N) = - G_final2(N+1:2*N, 1:N)/T;
% 
% minF = inf;
% bestalpha = 0;
% xx = 1:1:20;
% for ii = xx
%     hatalpha = 0.05 * ii;
%     Fn(ii) = norm(LLL_final2(1:N, 1:N) - LLL_final2(1:N, N+1:2*N)./hatalpha, 'fro');
%     if Fn(ii) < minF
%         minF = Fn(ii);
%         bestalpha = hatalpha;
%     end
% end
% L_calcu2 = (LLL_final2(1:N, 1:N) + LLL_final2(1:N, N+1:2*N)./bestalpha).*1/2;
% figure;
% plot(xx, Fn(1:end), 'LineWidth',2);hold on;grid on;


%% 绘图
%%% 静态 tvx 图
figure;
for i = 1:N
    plot(1:k_limits,global_x1(N + i,1:end-1),'LineWidth',2);xlim([0 k_limits]);ylim([-50 350]);hold on;grid on;
    set(gca,'FontSize',16);
    xlabel('Time-step','FontSize',16);ylabel('Velocities of the Agents','FontSize',16);
end
legend('Agent 1','Agent 2','Agent 3');
%%% 静态 tvy 图
figure;
for i = 1:N
    plot(1:k_limits,global_y1(N + i,1:end-1),'LineWidth',2);xlim([0 k_limits]);ylim([-250 250]);hold on;grid on;
    set(gca,'FontSize',16);
    xlabel('Time-step','FontSize',16);ylabel('Velocities of the Agents','FontSize',16);
end
legend('Agent 1','Agent 2','Agent 3');
%%% 静态 xy 图
figure;
for i = 1:N
    p1(i) = plot(global_x1(i,1:end-1) + targ_px(i) * ones(1, k_limits), global_y1(i,1:end-1) + targ_py(i) * ones(1, k_limits), 'LineWidth',2);
    xlim([0 3500]);ylim([800 2500]);hold on;grid on;
    set(gca,'FontSize',16);
    xlabel('X Axis','FontSize',16);ylabel('Y Axis','FontSize',16);
end
plot([global_x1(1,1) + targ_px(1), global_x1(2,1) + targ_px(2)], [global_y1(1,1) + targ_py(1), global_y1(2,1) + targ_py(2)], '-o', 'Color', 'b', 'LineWidth', 2);
plot([global_x1(3,1) + targ_px(3), global_x1(2,1) + targ_px(2)], [global_y1(3,1) + targ_py(3), global_y1(2,1) + targ_py(2)], '-o', 'Color', 'b', 'LineWidth', 2);
plot([global_x1(1,1) + targ_px(1), global_x1(3,1) + targ_px(3)], [global_y1(1,1) + targ_py(1), global_y1(3,1) + targ_py(3)], '-o', 'Color', 'b', 'LineWidth', 2);
plot([global_x1(1,end-1) + targ_px(1), global_x1(2,end-1) + targ_px(2)], [global_y1(1,end-1) + targ_py(1), global_y1(2,end-1) + targ_py(2)], '-o', 'Color', 'b', 'LineWidth', 2);
plot([global_x1(3,end-1) + targ_px(3), global_x1(2,end-1) + targ_px(2)], [global_y1(3,end-1) + targ_py(3), global_y1(2,end-1) + targ_py(2)], '-o', 'Color', 'b', 'LineWidth', 2);
plot([global_x1(1,end-1) + targ_px(1), global_x1(3,end-1) + targ_px(3)], [global_y1(1,end-1) + targ_py(1), global_y1(3,end-1) + targ_py(3)], '-o', 'Color', 'b', 'LineWidth', 2);
legend([p1(1),p1(2),p1(3)],'Agent 1','Agent 2','Agent 3','Location','Southeast');

%%% 动态 xy 图
% figure;
% color = [[0.00,0.45,0.74];[0.85,0.33,0.10];[0.93,0.69,0.13];[0.72,0.27,1];[0.47,0.67,0.19]];
% for k = 1:k_limits
%     for i = 1:N
%         plot([global_x1(i,k)+ targ_px(i),global_x1(i,k+1)+ targ_px(i)],[global_y1(i,k)+ targ_py(i),global_y1(i,k+1)+ targ_py(i)],'-.','MarkerSize',10,'Color',color(i,:));xlim([0 5000]);ylim([0 3000]);
%         hold on;
%     end
%     pause(0.05);   % 暂停，就可以看到点的变化走向
% end

%%% 静态 tvx 图
figure;
for i = 1:N
    plot(1:k_limits,global_x2(N + i,1:end-1),'LineWidth',2);xlim([0 k_limits]);ylim([-50 350]);hold on;grid on;
    set(gca,'FontSize',16);%xlim([0 k_limits]);ylim([0 350]);
    xlabel('Time-step','FontSize',16);ylabel('Velocities of the Agents','FontSize',16);
end
legend('Agent 1','Agent 2','Agent 3');

%%% 静态 tvy 图
figure;
for i = 1:N
    plot(1:k_limits,global_y2(N + i,1:end-1),'LineWidth',2);xlim([0 k_limits]);ylim([-250 250]);hold on;grid on;
    set(gca,'FontSize',16);% xlim([0 k_limits]);ylim([-250 250]);
    xlabel('Time-step','FontSize',16);ylabel('Velocities of the Agents','FontSize',16);
end
legend('Agent 1','Agent 2','Agent 3');

% %%% 静态 tv 图
% figure;
% for i = 1:N
%     plot(1:k_limits, sqrt(global_x2(N + i,1:end-1).^2 + global_y2(N + i,1:end-1).^2),'LineWidth',2);xlim([0 k_limits]);ylim([0 400]);hold on;grid on;
%     set(gca,'FontSize',16);
%     xlabel('Time-step','FontSize',16);ylabel('Velocities of the Agents','FontSize',16);
% end

%%% 静态 xy 图
figure;
for i = 1:N
    p2(i) = plot(global_x2(i,1:end-1) + targ_px(i) * ones(1, k_limits), global_y2(i,1:end-1) + targ_py(i) * ones(1, k_limits), 'LineWidth',2);
    xlim([0 3500]);ylim([800 2500]);hold on;grid on;
    set(gca,'FontSize',16);
    xlabel('X Axis','FontSize',16);ylabel('Y Axis','FontSize',16);
end
plot([global_x2(1,1) + targ_px(1), global_x2(2,1) + targ_px(2)], [global_y2(1,1) + targ_py(1), global_y2(2,1) + targ_py(2)], '-o', 'Color', 'b', 'LineWidth', 2);
plot([global_x2(3,1) + targ_px(3), global_x2(2,1) + targ_px(2)], [global_y2(3,1) + targ_py(3), global_y2(2,1) + targ_py(2)], '-o', 'Color', 'b', 'LineWidth', 2);
plot([global_x2(1,1) + targ_px(1), global_x2(3,1) + targ_px(3)], [global_y2(1,1) + targ_py(1), global_y2(3,1) + targ_py(3)], '-o', 'Color', 'b', 'LineWidth', 2);
plot([global_x2(1,end-1) + targ_px(1), global_x2(2,end-1) + targ_px(2)], [global_y2(1,end-1) + targ_py(1), global_y2(2,end-1) + targ_py(2)], '-o', 'Color', 'b', 'LineWidth', 2);
plot([global_x2(3,end-1) + targ_px(3), global_x2(2,end-1) + targ_px(2)], [global_y2(3,end-1) + targ_py(3), global_y2(2,end-1) + targ_py(2)], '-o', 'Color', 'b', 'LineWidth', 2);
plot([global_x2(1,end-1) + targ_px(1), global_x2(3,end-1) + targ_px(3)], [global_y2(1,end-1) + targ_py(1), global_y2(3,end-1) + targ_py(3)], '-o', 'Color', 'b', 'LineWidth', 2);
legend([p2(1),p2(2),p2(3)],'Agent 1','Agent 2','Agent 3','Location','Southeast');


%% 画图
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



