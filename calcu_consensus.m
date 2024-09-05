% 计算 k 时刻的 consensus 值
function [con_vx, con_vy, con_px, con_py] = calcu_consensus(L, init_vx, init_vy, init_px, init_py, targ_px, targ_py, k, T)
[n,~] = size(L);
trav_w = [zeros(1,n),ones(1,1)] * pinv([L,ones(n,1)]);
% 收敛速度
con_vx = trav_w * init_vx * ones(n, 1);
con_vy = trav_w * init_vy * ones(n, 1);
% 收敛位置（中间计算位置的时候用的是 relative position，所以这里要稍稍做处理）
% 另外 MATLAB 从 p(1) 开始，计算的是 p(k)，所以这里是 k-1
con_px = (trav_w * (init_px - targ_px) + (k - 1) * T * trav_w * init_vx) * ones(n, 1) + targ_px;
con_py = (trav_w * (init_py - targ_py) + (k - 1) * T * trav_w * init_vy) * ones(n, 1) + targ_py;
end