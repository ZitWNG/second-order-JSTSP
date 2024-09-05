function indicator = test_alphaT(L, alpha, T)
indicator = 0;
eigs = eig(-L);
[N,~] = size(L);
cnt = 0;
for i = 1:N
    eig_mu = eigs(i);
    real_mu = real(eig_mu);
    imag_mu = imag(eig_mu);
    abs_mu = abs(eig_mu);
    if real_mu < 0 && imag_mu ~= 0
        phi_mu = - 8 * imag_mu * imag_mu/abs_mu^4/(T-2*alpha)^2 - 2 * real_mu/abs_mu;
        if alpha * T < phi_mu
            cnt = cnt + 1;
        end
    end
    if real_mu < -0.0001 && imag_mu == 0
        phi_mu = -2/real_mu;
        if alpha * T < phi_mu
            cnt = cnt + 1;
        end
    end
end
if cnt == N-1
    indicator = 1;
end
end