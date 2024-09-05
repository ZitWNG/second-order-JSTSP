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
k_limits = 300;    % N^2;
repetition = 10;%10;

k_num = 300;
batchLength = 2;    % batchLength 条线做均值处理，处理掉图线上极端的 min max 值
batchNum = repetition/batchLength;
smoothwin = 1;


%% 普通迭代
dist_0 = zeros(repetition,k_num);
for r = 1:repetition
    [global_x0, global_y0] = inputDesign_ori(init_x, init_y, k_limits, G);
    for k = 1:k_num
        [~, G_final0] = estimation(global_x0(:,1:k_limits/k_num*k),G);
        L_calcu0 = (eye(N) - G_final0(1:N, 1:N)) * 2 / (T^2);
        dist_0(r,k) = norm(L_calcu0 - L, 2);
    end
end
summ0 = zeros(batchNum, k_num);
for i = 1:batchNum
%     summ0(i,:) = nansum(dist_0(batchLength * (i-1)+1:i*batchLength,:))/batchLength;
    summ0(i,:) = zeros(1,k_num);
end
y1 = smooth(max(summ0),smoothwin)';
y2 = smooth(min(summ0),smoothwin)';
y3 = smooth(sum(summ0)/batchNum,smoothwin)';


%% input design 1
seq1 = 1*[1 0 -3 2];
dist_1 = zeros(repetition,k_num);
for r = 1:repetition
    [global_x1, global_y1] = inputDesign_noisy(init_x, init_y, k_limits, G, T, seq1);
    for k = 1:k_num
        [~, G_final1] = estimation(global_x1, G);
        L_calcu1 = (eye(N) - G_final1(1:N, 1:N)) * 2 / (T^2);
        dist_1(r,k) = norm(L_calcu1 - L, 2);
    end
end
summ1 = zeros(batchNum, k_num);
for i = 1:batchNum
    summ1(i,:) = sum(dist_1(batchLength * (i-1)+1:i*batchLength,:))/batchLength;
end
y4 = smooth(max(summ1),smoothwin)';
y5 = smooth(min(summ1),smoothwin)';
y6 = smooth(sum(summ1)/batchNum,smoothwin)';

%% input design 2
seq2 = 10*[1 0 -3 2];
dist_2 = zeros(repetition,k_num);
for r = 1:repetition
    [global_x2, global_y2] = inputDesign_noisy(init_x, init_y, k_limits, G, T, seq2);
    for k = 1:k_num
        [~, G_final2] = estimation(global_x2, G);
        L_calcu2 = (eye(N) - G_final2(1:N, 1:N)) * 2 / (T^2);
        dist_2(r,k) = norm(L_calcu2 - L, 2);
    end
end
summ2 = zeros(batchNum, k_num);
for i = 1:batchNum
    summ2(i,:) = sum(dist_2(batchLength * (i-1)+1:i*batchLength,:))/batchLength;
end
y7 = smooth(max(summ2),smoothwin)';
y8 = smooth(min(summ2),smoothwin)';
y9 = smooth(sum(summ2)/batchNum,smoothwin)';

%% input design 3
seq3 = 10*[1 -1.5 0 0.5];
dist_3 = zeros(repetition,k_num);
for r = 1:repetition
    [global_x3, global_y3] = inputDesign_noisy(init_x, init_y, k_limits, G, T, seq3);
    for k = 1:k_num
        [~, G_final3] = estimation(global_x3, G);
        L_calcu3 = (eye(N) - G_final3(1:N, 1:N)) * 2 / (T^2);
        dist_3(r,k) = norm(L_calcu3 - L, 2);
    end
end
summ3 = zeros(batchNum, k_num);
for i = 1:batchNum
    summ3(i,:) = sum(dist_3(batchLength * (i-1)+1:i*batchLength,:))/batchLength;
end
y10 = smooth(max(summ3),smoothwin)';
y11 = smooth(min(summ3),smoothwin)';
y12 = smooth(sum(summ3)/batchNum,smoothwin)';

%% figure
figure;
ccolor = [[152,78,163];[55,126,184];[77,175,74];[228,26,28]]./256; % 紫蓝绿红
klist = k_limits/k_num:k_limits/k_num:k_limits;
area1 = fill([klist fliplr(klist)],[y1 fliplr(y2)],ccolor(1,:),'edgealpha', '0', 'facealpha', '.5');hold on
area2 = fill([klist fliplr(klist)],[y4 fliplr(y5)],ccolor(2,:),'edgealpha', '0', 'facealpha', '.5');hold on
area3 = fill([klist fliplr(klist)],[y7 fliplr(y8)],ccolor(3,:),'edgealpha', '0', 'facealpha', '.5');hold on
area4 = fill([klist fliplr(klist)],[y10 fliplr(y11)],ccolor(4,:),'edgealpha', '0', 'facealpha', '.5');hold on
plot(klist,y3,'Color',ccolor(1,:),'LineWidth',1.5);hold on
plot(klist,y6,'Color',ccolor(2,:),'LineWidth',1.5);hold on
plot(klist,y9,'Color',ccolor(3,:),'LineWidth',1.5);hold on
plot(klist,y12,'Color',ccolor(4,:),'LineWidth',1.5);
% values = spcrv([[k_limits_list(1) k_limits_list k_limits_list(end)];[y6(1) y6 y6(end)]],10);
% plot(values(1,:),values(2,:), 'b');

% lgd = legend('Error Range using normal protocol','Error Range using \{1,0,-3,2\}','Error Range using \{10,0,-30,20\}','Error Range using \{10,-15,0,5\}','Average Error using normal protocol','Average Error using \{10,0,-30,20\}','Average Error using \{1,0,-3,2\}','Average Error using \{10,-15,0,5\}');
lgd = legend([area1, area2, area3, area4],'Error Range using normal protocol','Error Range using \{1,0,-3,2\}','Error Range using \{10,0,-30,20\}','Error Range using \{10,-15,0,5\}','Location','North');
lgd.FontSize = 16;
set(gca,'FontSize',18);
xlabel('Observation Time k','FontSize',24);ylabel('Inference Error','FontSize',24);
axis([0 300 -5 200]);

%%



