function SeMin=objfun(para,x)
% para(1)-alpha �ͷ�����
% para(2)-K ��������
alpha=round(para(1));
K=round(para(2));

tau = 0;            % noise-tolerance (no strict fidelity enforcement)
DC = 0;             % no DC part imposed
init = 1;           % initialize omegas uniformly
tol = 1e-5;
%% VMD�ֽ�
[u, u_hat, omega] = VMD(x, alpha, tau, K, DC, init, tol);
%% ��ȡ������
[m,n]=size(u);
mm=2;
for ii=1:m
feature(ii)=SampEn(u(ii,:), mm, 0.2*std(u(ii,:)));%u�������أ�0.2*std(imf1(1,:))��ʾ���������r��ֵ��
end
SeMin=min(feature);