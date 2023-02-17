function SeMin=objfun(para,x)
% para(1)-alpha 惩罚参数
% para(2)-K 分量个数
alpha=round(para(1));
K=round(para(2));

tau = 0;            % noise-tolerance (no strict fidelity enforcement)
DC = 0;             % no DC part imposed
init = 1;           % initialize omegas uniformly
tol = 1e-5;
%% VMD分解
[u, u_hat, omega] = VMD(x, alpha, tau, K, DC, init, tol);
%% 提取样本熵
[m,n]=size(u);
mm=2;
for ii=1:m
feature(ii)=SampEn(u(ii,:), mm, 0.2*std(u(ii,:)));%u的样本熵，0.2*std(imf1(1,:))表示求解样本熵r阀值，
end
SeMin=min(feature);