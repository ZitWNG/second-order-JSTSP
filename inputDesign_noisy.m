%%% noisy input second-order consensus algorithm
function [global_x2, global_y2] = inputDesign_noisy(init_x, init_y, k_limits, G, T, seq)
global_x2 = init_x;
global_y2 = init_y;
[n,~] = size(G);
N = n / 2;
b = rand(N, k_limits);      % start indicator
% seq = [1 -1 -1 1];  % [1 -2 1];
[~, lens] = size(seq);
theta_x2 = zeros(N, k_limits + lens);
theta_y2 = zeros(N, k_limits + lens);
a_x2 = zeros(N, k_limits);
a_y2 = zeros(N, k_limits);

amp = 100;
phi = 0.99;

for k = 1:k_limits
    for i = 1:N
        if b(i, k) > 0.5 % k < k_limits - 40 && b(i, k) > 0.5
            a_x2(i, k) = amp * phi^k * rand(1);
            a_y2(i, k) = amp * phi^k * rand(1);
            for ell = 1:lens
                theta_x2(i, k + ell - 1) = seq(1, ell) * a_x2(i, k) + theta_x2(i, k + ell - 1);
                theta_y2(i, k + ell - 1) = seq(1, ell) * a_y2(i, k) + theta_y2(i, k + ell - 1);
            end
        end
    end
    input_x2 = [T^2/2 * theta_x2(:, k); T * theta_x2(:, k)];
    input_y2 = [T^2/2 * theta_y2(:, k); T * theta_y2(:, k)];
    global_x2(:,k+1) = G * global_x2(:,k) + input_x2;
    global_y2(:,k+1) = G * global_y2(:,k) + input_y2;
end
end

