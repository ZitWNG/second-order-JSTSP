%%% no input
function [global_x] = inputDesign_0(init_state, k_limits, G, T)
global_x = init_state;
[n,~] = size(G);
N = n / 2;
for k = 1:k_limits
    global_x(:,k+1) = G * global_x(:,k);
end

end