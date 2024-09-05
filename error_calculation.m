%% 计算估计的P矩阵与真实P矩阵误差
% P实际的P，P_est估计的P，B判断连接矩阵，B_est估计的B，flag取不同数值对应不同的误差求解
function [error_sum] = error_calculation(P,P_est,flag)
[N,~] = size(P);

if (flag == 2)                  % 结构误差(连接与否)
    B = sign(P);
    P(P == 0) = inf;
    bound = 0.5 * min(min(P));  % why 0.5
    P_est(P_est <= bound) = 0;  % 小于阈值的元素置零（小于阈值置零是合理的，但是P_es并没有该先验知识，而且别的方法都没有进行类似的处理）
    B_est = sign(P_est);        % 二值化得到连接矩阵
    error_sum = sum(sum(abs(B_est-B)));
    error_sum = error_sum/(N^2);
elseif (flag == 1)              % 2-norm 误差
    error_sum = norm(P - P_est);
    error_sum = error_sum/norm(P);
elseif(flag == 0)               % 幅值误差(F-norm)
    error_sum = norm(P-P_est, 'fro' );
    error_sum = error_sum/norm(P, 'fro' );
end

end