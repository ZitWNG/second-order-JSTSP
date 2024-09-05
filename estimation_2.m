function [output, L_final] = estimation_2(states, L_real, T, alpha)
%%% states为一致性过程运行的状态轨迹  W_real为交互矩阵

%% main parameters
%%%% regression method
method_flag = 1;%4;    % 1求逆 2 Lasso  3 NLLS  4 CLLS
lambda = 3;         % Lasso regularization parameter

[N,~] = size(states);
n = N / 2;
%%% 额外构造 X Y 数据
x_1a = states(1:n, 1:end-1);
x_1b = states(1:n, 2:end);
x_2a = states(n+1:N, 1:end-1);
x_2b = states(n+1:N, 2:end);
Za_1 = T^2/2 * x_1a + alpha * T^2/2 * x_2a;
Za_2 = T * x_1a + alpha * T * x_2a;
Zb_1 = x_1a + T * x_2a - x_1b;
Zb_2 = x_2a - x_2b;
Za = [Za_1, Za_2];
Zb = [Zb_1, Zb_2];

%% 转化为向量 带求解的未知向量是 W_temp (也是按行顺序拼接排列的一列向量)
%%% 构造对角块矩阵Za_temp
Za_temp = Za';
for j = 1:n-1
    Za_temp = blkdiag(Za_temp,Za');
end
%%% 构造块向量Zb_temp
Zb_temp = Zb';
Zb_temp = Zb_temp(:); % 直接得到 Zb_temp 矩阵按行顺序拼接排列的一列向量
L_temp = get_L(Za_temp, Zb_temp, n, lambda, method_flag);
L_final = reshape(L_temp, n, n);
L_final = L_final';

error_s_sum = [];
error_s_sum(1,1) = error_calculation(L_real,L_final,1);    % 2-norm 误差
error_s_sum(1,2) = error_calculation(L_real,L_final,0);    % 幅值误差(F-norm)

output = [error_s_sum(1,1),error_s_sum(1,2)];

end