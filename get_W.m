% 求矩阵
% Inputs:   Za_temp - 系数矩阵
%           Zb_temp - measurements with noise
%           N - number ofagents
%           lambda - Lasso regularization parameter
%           method_flag - 求解方法
% Outputs: P_final - 矩阵
function [P_final] = get_W(Za_temp,Zb_temp,N,lambda,method_flag)
C = Za_temp;
y = Zb_temp;  % measurements with noise
if(method_flag == 1)
%% 方案一 最小二乘 拼接成一列向量1次求逆
    P_temp = inv(Za_temp'*Za_temp)*Za_temp'*Zb_temp;
elseif(method_flag == 2)
%% 方案二 用Lasso优化方法直接求 有非负约束和正则项, min|Za_temp P_temp-Zb_temp|^2+lamda|P_temp|_1
    %%%%优化参数 lambda 的选取非常关键！1~10之间比较合适
    % lambda = 3;      % regularization parameter
    rel_tol = 0.001;     % relative target duality gap
    [P_temp,~] = l1_ls_nonneg(C,y,lambda,rel_tol,1);
elseif(method_flag == 3)
%% 方案三 matlab自带求解器 Nonnegative linear-least-squares(NLLS)
    P_temp = lsqnonneg(C,y);
elseif(method_flag == 4)
%% 方案四 matlab自带求解器 Constrained linear-least-squares(CLLS)
%%%%% 构造单位约束矩阵及系数 P1=1  如果P是向量 则1P=1
    U1 = ones(1,N);
    for i = 1:N-1
        U1 = blkdiag(U1,ones(1,N));
    end
    U2 = ones(N,1); %%%直接得到Zb_temp矩阵按行顺序拼接排列的一列向
    [num,~] = size(y);
    lb = zeros(num,1); %%下界为0
    ub = ones(num,1);  %%上界为1
    P_temp = lsqlin(C,y,[],[],U1,U2,lb,ub); %%% 标准用法lsqlin(C,d,A,b,Aeq,beq,lb,ub) CLLS  P拼接成一列算1次
end
P_final = P_temp;
%% 将得到的向量P_temp转化为对应的矩阵
%%% reshape函数是按列排列的, 第一个参数N表示保留的行数，第二个参数N表示保留的列数 
% P_final=reshape(P_temp,N,N); 
% P_final=P_final';
end
