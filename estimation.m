function [output, W_final] = estimation(states,W_real)
%%% states为一致性过程运行的状态轨迹  W_real为交互矩阵

%% main parameters
%%%% regression method
method_flag = 1;%4;    % 1求逆 2 Lasso  3 NLLS  4 CLLS
lambda = 3;         % Lasso regularization parameter

[N,~] = size(states);

Za = states(:,1:end-1);
Zb = states(:,2:end);
%% 转化为向量 带求解的未知向量是 W_temp (也是按行顺序拼接排列的一列向量)
%%% 构造对角块矩阵Za_temp
Za_temp = Za';
for j = 1:N-1
    Za_temp = blkdiag(Za_temp,Za');
end
%%% 构造块向量Zb_temp
Zb_temp = Zb';
Zb_temp = Zb_temp(:); %%%直接得到Zb_temp矩阵按行顺序拼接排列的一列向量
W_temp = get_W(Za_temp,Zb_temp,N,lambda,method_flag);
W_final = reshape(W_temp,N,N);
W_final = W_final';

error_s_sum = [];
error_s_sum(1,1) = error_calculation(W_real,W_final,1);    % 2-norm 误差
error_s_sum(1,2) = error_calculation(W_real,W_final,0);    % 幅值误差(F-norm)

output = [error_s_sum(1,1),error_s_sum(1,2)];

end