% �����
% Inputs:   Za_temp - ϵ������
%           Zb_temp - measurements with noise
%           N - number ofagents
%           lambda - Lasso regularization parameter
%           method_flag - ��ⷽ��
% Outputs: P_final - ����
function [P_final] = get_W(Za_temp,Zb_temp,N,lambda,method_flag)
C = Za_temp;
y = Zb_temp;  % measurements with noise
if(method_flag == 1)
%% ����һ ��С���� ƴ�ӳ�һ������1������
    P_temp = inv(Za_temp'*Za_temp)*Za_temp'*Zb_temp;
elseif(method_flag == 2)
%% ������ ��Lasso�Ż�����ֱ���� �зǸ�Լ����������, min|Za_temp P_temp-Zb_temp|^2+lamda|P_temp|_1
    %%%%�Ż����� lambda ��ѡȡ�ǳ��ؼ���1~10֮��ȽϺ���
    % lambda = 3;      % regularization parameter
    rel_tol = 0.001;     % relative target duality gap
    [P_temp,~] = l1_ls_nonneg(C,y,lambda,rel_tol,1);
elseif(method_flag == 3)
%% ������ matlab�Դ������ Nonnegative linear-least-squares(NLLS)
    P_temp = lsqnonneg(C,y);
elseif(method_flag == 4)
%% ������ matlab�Դ������ Constrained linear-least-squares(CLLS)
%%%%% ���쵥λԼ������ϵ�� P1=1  ���P������ ��1P=1
    U1 = ones(1,N);
    for i = 1:N-1
        U1 = blkdiag(U1,ones(1,N));
    end
    U2 = ones(N,1); %%%ֱ�ӵõ�Zb_temp������˳��ƴ�����е�һ����
    [num,~] = size(y);
    lb = zeros(num,1); %%�½�Ϊ0
    ub = ones(num,1);  %%�Ͻ�Ϊ1
    P_temp = lsqlin(C,y,[],[],U1,U2,lb,ub); %%% ��׼�÷�lsqlin(C,d,A,b,Aeq,beq,lb,ub) CLLS  Pƴ�ӳ�һ����1��
end
P_final = P_temp;
%% ���õ�������P_tempת��Ϊ��Ӧ�ľ���
%%% reshape�����ǰ������е�, ��һ������N��ʾ�������������ڶ�������N��ʾ���������� 
% P_final=reshape(P_temp,N,N); 
% P_final=P_final';
end
