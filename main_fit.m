%% MYT 和 YYT-1
n1 = zeros(1,300);
n2 = zeros(1,300);
n3 = zeros(1,300);
m1 = zeros(1,300);
m2 = zeros(1,300);
m3 = zeros(1,300);
for x=10:10:3000
%     y(1,x/100) = norm(Theta_2(:,1:x) * Y2(:,1:x)'*inv(Y2(:,1:x)*Y2(:,1:x)'), 'fro');
    n1(1,x/10) = norm(inv(Y1(:,1:x)*Y1(:,1:x)'), 'fro');
    n2(1,x/10) = norm(inv(Y2(:,1:x)*Y2(:,1:x)'), 'fro');
    n3(1,x/10) = norm(inv(Y3(:,1:x)*Y3(:,1:x)'), 'fro');
    m1(1,x/10) = norm(Theta_1(:,1:x) * Y1(:,1:x)', 'fro');
    m2(1,x/10) = norm(Theta_2(:,1:x) * Y2(:,1:x)', 'fro');
    m3(1,x/10) = norm(Theta_3(:,1:x) * Y3(:,1:x)', 'fro');
end
x = 10:10:3000;
figure;
plot(x,n1,'b');hold on
plot(x,n2,'g');hold on
plot(x,n3,'r');
figure;
plot(x,m1,'b');hold on
plot(x,m2,'g');hold on
plot(x,m3,'r');


%% 拟合图线
% modelfun = @(b,x)(b(1)+b(2)*exp(-b(3)*x));
% modelfun = @(b,x)(1./(b(1) + b(2)*exp(-b(3)*x) + b(4)*x));
modelfun = @(b,x)(1./(b(1) + b(2)*x));
b0 = [0.001 0.001];
beta = nlinfit(x, n1, modelfun, b0);
figure;
plot(x,modelfun(beta,x))

modelfun = @(b,x)(b(1) + b(2) * x);
b0 = [0.001 0.001];
beta = nlinfit(x, m1, modelfun, b0);
figure;
plot(x,modelfun(beta,x))

modelfun = @(b,x)((b(1) + b(2) * x)./(b(3) + b(4) * x));
b0 = [0.001 0.001 0.001 0.001];
yy1 = n1.*m1;
beta = nlinfit(x, yy1, modelfun, b0);
figure;
plot(x,modelfun(beta,x))


%%
