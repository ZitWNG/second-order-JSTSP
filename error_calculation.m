%% ������Ƶ�P��������ʵP�������
% Pʵ�ʵ�P��P_est���Ƶ�P��B�ж����Ӿ���B_est���Ƶ�B��flagȡ��ͬ��ֵ��Ӧ��ͬ��������
function [error_sum] = error_calculation(P,P_est,flag)
[N,~] = size(P);

if (flag == 2)                  % �ṹ���(�������)
    B = sign(P);
    P(P == 0) = inf;
    bound = 0.5 * min(min(P));  % why 0.5
    P_est(P_est <= bound) = 0;  % С����ֵ��Ԫ�����㣨С����ֵ�����Ǻ���ģ�����P_es��û�и�����֪ʶ�����ұ�ķ�����û�н������ƵĴ���
    B_est = sign(P_est);        % ��ֵ���õ����Ӿ���
    error_sum = sum(sum(abs(B_est-B)));
    error_sum = error_sum/(N^2);
elseif (flag == 1)              % 2-norm ���
    error_sum = norm(P - P_est);
    error_sum = error_sum/norm(P);
elseif(flag == 0)               % ��ֵ���(F-norm)
    error_sum = norm(P-P_est, 'fro' );
    error_sum = error_sum/norm(P, 'fro' );
end

end