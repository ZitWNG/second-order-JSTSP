% �Աȼ�������� ������ÿ��ʱ�̵�״̬����״̬ƽ����Ļع�Ч���Ա�
function [output] = regression_with_noise(states,W_real)
noise_v = 0:0.1:0.5;
% noise_v = 0:0.2:1;
noise_v = noise_v';
[noise_num,~] = size(noise_v);

%% �ؼ������趨
%%%%�ع鷽��ѡ��
% method_flag = 3;  %%%%% ѡ��������ⷽ�� =1 ���� 2 Lasso  3 NLLS  4 CLLS
% lambda=3; %%lambda��Lasso regularization parameter    
random_num = 1;  %%% ͬһ����ǿ�������ʵ�����

[N,countmax]=size(states);

for q = 1:noise_num
    %% ��ԭʼ���ݼӹ۲�����
    sigma = noise_v(q); % ����ǿ��
    %% ���趨����ģʽ���random_num�μ�����ƽ��Ч��
    error_sum = []; % �ϴγ������к����û����ñ��� Ҫ��ʼ�����

    %% ÿ������ǿ���¼�random_num�κ����ƽ��
    for i = 1:random_num
        z_states = states + normrnd(0,sigma,N,countmax);
        error_sum(i,:) = get_interaction(z_states, W_real); %% ��i��ͬһ����ǿ���µ�������������
        
        %%% ƽ���ķ���
        [new_states] = ave_state( z_states );
         new_error_sum(i,:)= get_interaction(new_states,W_real); %% ��i��ͬһ����ǿ���µ�������������
        
    end
    %% random_num�ε�ͬһ����ǿ���µĻع�ƽ�����
    uu = mean(error_sum,1);%%ÿ����Ps�����
    er_layer11(1,q) = uu(1,1);  %%ֻȡ��1һ��
    er_layer12(1,q) = uu(1,2);  %%ֻȡ��1һ��
    
    %% random_num�ε�ͬһ����ǿ���µĻع�ƽ����� 
    new_uu = mean(new_error_sum,1);%%ÿ����Ps�����
    new_er_layer11(1,q) = new_uu(1,1);  %%ֻȡ��1һ��
    new_er_layer12(1,q) = new_uu(1,2);  %%ֻȡ��1һ��
    
end

% estimated_

output=[er_layer11;er_layer12];

%% ��ͼ
figure;
%%% ����Ȩ
plot(1:noise_num,er_layer11,'-^b','LineWidth',1.5);hold on;
plot(1:noise_num,er_layer12,'-sr','LineWidth',1.5);hold on;

%%% ��״̬ƽ���Ľ��
plot(1:noise_num,new_er_layer11,'-.^b','LineWidth',1.5);hold on;
plot(1:noise_num,new_er_layer12,'-.sr','LineWidth',1.5);hold on;

grid on;
%%%% x����ʾ��Ӧ�ľ�������ֵ ע��num2str(noise_v)��noise_v����Ϊһ������
xticks(1:1:noise_num);
xticklabels( num2str(noise_v) );
xlabel('Noise Variance','FontSize',14,'Vertical','top');
ylabel('Recovery Errors','FontSize',14,'Vertical','middle');
title('The approximation error of the internal interaction rule','FontSize',14);
h1=legend('Structure error','Magnitude error','Structure error (a)','Magnitude error(b)');
set(h1,'FontSize',14);
set(gca,'FontSize',14);

end

