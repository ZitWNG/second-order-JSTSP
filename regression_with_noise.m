% 对比加噪情况下 纯依赖每个时刻的状态与多个状态平均后的回归效果对比
function [output] = regression_with_noise(states,W_real)
noise_v = 0:0.1:0.5;
% noise_v = 0:0.2:1;
noise_v = noise_v';
[noise_num,~] = size(noise_v);

%% 关键参数设定
%%%%回归方法选择
% method_flag = 3;  %%%%% 选择哪种求解方法 =1 求逆 2 Lasso  3 NLLS  4 CLLS
% lambda=3; %%lambda是Lasso regularization parameter    
random_num = 1;  %%% 同一噪声强度下随机实验次数

[N,countmax]=size(states);

for q = 1:noise_num
    %% 给原始数据加观测噪声
    sigma = noise_v(q); % 噪声强度
    %% 按设定参数模式求解random_num次加噪后的平均效果
    error_sum = []; % 上次程序运行后可能没清楚该变量 要初始化清空

    %% 每个噪声强度下加random_num次后误差平均
    for i = 1:random_num
        z_states = states + normrnd(0,sigma,N,countmax);
        error_sum(i,:) = get_interaction(z_states, W_real); %% 第i次同一噪声强度下的随机加噪误差结果
        
        %%% 平均的方法
        [new_states] = ave_state( z_states );
         new_error_sum(i,:)= get_interaction(new_states,W_real); %% 第i次同一噪声强度下的随机加噪误差结果
        
    end
    %% random_num次的同一噪声强度下的回归平均误差
    uu = mean(error_sum,1);%%每列是Ps的误差
    er_layer11(1,q) = uu(1,1);  %%只取第1一列
    er_layer12(1,q) = uu(1,2);  %%只取第1一列
    
    %% random_num次的同一噪声强度下的回归平均误差 
    new_uu = mean(new_error_sum,1);%%每列是Ps的误差
    new_er_layer11(1,q) = new_uu(1,1);  %%只取第1一列
    new_er_layer12(1,q) = new_uu(1,2);  %%只取第1一列
    
end

% estimated_

output=[er_layer11;er_layer12];

%% 画图
figure;
%%% 不加权
plot(1:noise_num,er_layer11,'-^b','LineWidth',1.5);hold on;
plot(1:noise_num,er_layer12,'-sr','LineWidth',1.5);hold on;

%%% 带状态平均的结果
plot(1:noise_num,new_er_layer11,'-.^b','LineWidth',1.5);hold on;
plot(1:noise_num,new_er_layer12,'-.sr','LineWidth',1.5);hold on;

grid on;
%%%% x轴显示对应的具体噪声值 注意num2str(noise_v)中noise_v必须为一列数组
xticks(1:1:noise_num);
xticklabels( num2str(noise_v) );
xlabel('Noise Variance','FontSize',14,'Vertical','top');
ylabel('Recovery Errors','FontSize',14,'Vertical','middle');
title('The approximation error of the internal interaction rule','FontSize',14);
h1=legend('Structure error','Magnitude error','Structure error (a)','Magnitude error(b)');
set(h1,'FontSize',14);
set(gca,'FontSize',14);

end

