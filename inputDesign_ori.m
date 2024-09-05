%%% original second-order consensus algorithm
function [global_x1, global_y1] = inputDesign_ori(init_x, init_y, k_limits, G)
global_x1 = init_x;
global_y1 = init_y;

for k = 1:k_limits
    global_x1(:,k+1) = G * global_x1(:,k);
    global_y1(:,k+1) = G * global_y1(:,k);
end

end
