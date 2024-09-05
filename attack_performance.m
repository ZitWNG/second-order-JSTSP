function [distortion, deviation] = attack_performance(con_vx, con_vy, con_px, con_py, global_x, global_y, targ_px, targ_py, alpha)
[n,~] = size(global_x);
N = n / 2;
last_px = global_x(1:N,end);last_py = global_y(1:N,end);                            % 最后时刻的 relative position
last_vx = global_x(N+1:end,end);last_vy = global_y(N+1:end,end);                    % 最后时刻的 velocity
delta_p = sqrt((con_px - targ_px - last_px).^2 + (con_py - targ_py - last_py).^2);  % 与 consensus position 的差异
delta_v = sqrt((con_vx - last_vx).^2 + (con_vy - last_vy).^2);                      % 与 consensus velocity 的差异
deviation = sum(delta_p) + alpha * sum(delta_v);
distortion_v = sqrt((last_px - sum(last_px)/N * ones(N,1)).^2 + (last_py - sum(last_py)/N * ones(N,1)).^2);
distortion = sum(abs(distortion_v));
end