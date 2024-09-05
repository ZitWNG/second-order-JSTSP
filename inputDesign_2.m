%%% SE 2
function [global_x, theta] = inputDesign_2(init_state, k_limits, G, T, phi)
global_x = init_state;
[n,~] = size(G);
N = n/2;
b = rand(N, k_limits);      % start indicator
theta = zeros(N, k_limits);
for k = 1:k_limits
    for i = 1:N
        if k < k_limits - 10 && b(i, k) > 0.1 %k < k_limits / 2 && b(i, k) > 0.9
            theta(i, k) = phi^k + theta(i, k);
            Tm = 1;%randi([1 5]);  % compensation period
            Tn = 3;%randi([Tm + 1 Tm + 5]);
            theta(i, k + Tm) = -Tn / (Tn-Tm) * phi^k + theta(i, k + Tm);
            theta(i, k + Tn) = Tm / (Tn-Tm) * phi^k + theta(i, k + Tn);
        end
    end
    input = [T^2/2 * theta(:, k); T * theta(:, k)];
    global_x(:,k+1) = G * global_x(:,k) + input;
end

end