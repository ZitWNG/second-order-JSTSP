%%% SE 1
function [global_x, theta] = inputDesign_1(init_state, k_limits, G, T, phi, Tc)
global_x = init_state;
[n,~] = size(G);
N = n/2;
b = rand(N, k_limits);      % start indicator 
theta = zeros(N, k_limits);
for k = 1:k_limits
    for i = 1:N
        if k < k_limits - 10 && b(i, k) > 0.5
            theta(i, k) = phi^k + theta(i, k);
            % Tc = randi([2 5]);          % compensation period
            theta(i, k + Tc) = phi^k + theta(i, k + Tc);
            for l = 1:Tc - 1
                theta(i, k + l) = -2 * phi^k / (Tc - 1) + theta(i, k + l);
            end
        end
    end
    input = [T^2/2 * theta(:, k); T * theta(:, k)];
    global_x(:,k+1) = G * global_x(:,k) + input;
end

end