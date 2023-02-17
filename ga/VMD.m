function [u, u_hat, omega] = VMD(signal, alpha, tau, K, DC, init, tol)
% Variational Mode Decomposition VMD变分模态分解  经验模式分解（EMD）
%EMD分解经常被用作信号特征提取的一个预先处理手段，将各IMF分量作为后续分析方法的输入，以完成更加复杂的工作。
% Authors: Konstantin Dragomiretskiy and Dominique Zosso
% zosso@math.ucla.edu --- http://www.math.ucla.edu/~zosso
% Initial release 2013-12-12 (c) 2013
%
% Input and Parameters输入和参数:
% ---------------------
% signal  - 此参数为一个信号，从表中读取的参数 the time domain signal (1D) to be decomposed
% -----alpha   -输入为     带宽限制，经验取值为抽样点长度1.5-2.0 倍；数据保真约束的平衡参数the balancing parameter of the data-fidelity constraint
% tau     - 输入为0    双爬升的时间步长(噪声松弛选0)time-step of the dual ascent ( pick 0 for noise-slack )
% -----K       - 输入为       需要恢复的模式个数，分解模态（IMF）个数the number of modes to be recovered
%DC      - 输入为0   如果将第一模态设置为DC(0频率)则为真true if the first mode is put and kept at DC (0-freq)
% init    - 输入为1      0 = all omegas start at 0  从零开始
%                    1 = all omegas start uniformly distributed   均匀分布
%                    2 = all omegas initialized randomly  随机分布
% tol     -输入为1e-5     收敛准则的容差;通常大约e-6tolerance of convergence criterion; typically around 1e-6
%迭代次数是反复迭代的次数，迭代次数越大，结果越准确。容差是计算结果和设定结果（期望结果）之间的误差值。

% Output:
% -------
% u       -分解模式的集合 the collection of decomposed modes
% u_hat   - 模谱spectra of the modes
% omega   - 估计模式中心频率 estimated mode center-frequencies
%
% When using this code, please do cite our paper:
% -----------------------------------------------
% K. Dragomiretskiy, D. Zosso, Variational Mode Decomposition, IEEE Trans.
% on Signal Processing (in press)
% please check here for update reference: 
%          http://dx.doi.org/10.1109/TSP.2013.2288675



%---------- Preparations
%signal=load('moni_noise.dat');alpha=3.216578811729743e+03;tau=0;K=5;DC=0;init=1;tol=1e-5

% 输入信号的周期和采样频率Period and sampling frequency of input signal
save_T = length(signal);%信号的时间2048
fs = 1/save_T;%信号的频率

% 通过镜像扩展信号extend the signal by mirroring
T = save_T;
f_mirror(1:T/2) = signal(T/2:-1:1);
f_mirror(T/2+1:3*T/2) = signal;
f_mirror(3*T/2+1:2*T) = signal(T:-1:T/2+1);
f = f_mirror;

% 时域0到T(镜像信号)Time Domain 0 to T (of mirrored signal)
T = length(f);%信号的长度4096，扩大了一倍
t = (1:T)/T;

% Spectral Domain discretization
freqs = t-0.5-1/T;

% 最大的迭代次数(如果还没有收敛，那么它也不会收敛)。Maximum number of iterations (if not converged yet, then it won't anyway)
N = 500;

% 为了将来的推广:每个模式都有单独的alphaFor future generalizations: individual alpha for each mode
Alpha = alpha*ones(1,K);

% Construct and center f_hat
f_hat = fftshift((fft(f)));
f_hat_plus = f_hat;
f_hat_plus(1:T/2) = 0;

% matrix keeping track of every iterant // could be discarded for mem
u_hat_plus = zeros(N, length(freqs), K);

%omega_k的初始化 Initialization of omega_k
omega_plus = zeros(N, K);
switch init%init输入为1
    case 1
        for i = 1:K
            omega_plus(1,i) = (0.5/K)*(i-1);
        end
    case 2
        omega_plus(1,:) = sort(exp(log(fs) + (log(0.5)-log(fs))*rand(1,K)));
    otherwise
        omega_plus(1,:) = 0;
end

% 如果施加直流模式，将其ω设为0if DC mode imposed, set its omega to 0
if DC%DC==0
    omega_plus(1,1) = 0;
end

% 从空的二元变量开始。start with empty dual variables
lambda_hat = zeros(N, length(freqs));

% other inits
uDiff = tol+eps; % update step
n = 1; %循环计数器 loop counter
sum_uk = 0; % accumulator



% ----------- Main loop for iterative updates主循环用于迭代更新
while ( uDiff > tol &&  n < N ) % 未收敛且低于迭代极限not converged and below iterations limit
    
    % update first mode accumulator
    k = 1;
    sum_uk = u_hat_plus(n,:,K) + sum_uk - u_hat_plus(n,:,1);
    
    % update spectrum of first mode through Wiener filter of residuals
    u_hat_plus(n+1,:,k) = (f_hat_plus - sum_uk - lambda_hat(n,:)/2)./(1+Alpha(1,k)*(freqs - omega_plus(n,k)).^2);
    
    % update first omega if not held at 0
    if ~DC%~DC==1
        omega_plus(n+1,k) = (freqs(T/2+1:T)*(abs(u_hat_plus(n+1, T/2+1:T, k)).^2)')/sum(abs(u_hat_plus(n+1,T/2+1:T,k)).^2);
    end
    
    % update of any other mode
    for k=2:K
        
        % accumulator
        sum_uk = u_hat_plus(n+1,:,k-1) + sum_uk - u_hat_plus(n,:,k);
        
        % 频率mode spectrum
        u_hat_plus(n+1,:,k) = (f_hat_plus - sum_uk - lambda_hat(n,:)/2)./(1+Alpha(1,k)*(freqs - omega_plus(n,k)).^2);
        
        %中心频率 center frequencies
        omega_plus(n+1,k) = (freqs(T/2+1:T)*(abs(u_hat_plus(n+1, T/2+1:T, k)).^2)')/sum(abs(u_hat_plus(n+1,T/2+1:T,k)).^2);
        
    end
    
    % Dual ascent
    lambda_hat(n+1,:) = lambda_hat(n,:) + tau*(sum(u_hat_plus(n+1,:,:),3) - f_hat_plus);
    
    % loop counter
    n = n+1;
    
    % converged yet?
    uDiff = eps;
    for i=1:K
        uDiff = uDiff + 1/T*(u_hat_plus(n,:,i)-u_hat_plus(n-1,:,i))*conj((u_hat_plus(n,:,i)-u_hat_plus(n-1,:,i)))';
    end
    uDiff = abs(uDiff);
    
end


%------ 后处理和清理Postprocessing and cleanup


% 如果提前聚合，则丢弃空空间discard empty space if converged early
N = min(N,n);
omega = omega_plus(1:N,:);%最终输出变量    估计模式中心频率

%信号重建 Signal reconstruction
u_hat = zeros(T, K);
u_hat((T/2+1):T,:) = squeeze(u_hat_plus(N,(T/2+1):T,:));
u_hat((T/2+1):-1:2,:) = squeeze(conj(u_hat_plus(N,(T/2+1):T,:)));
u_hat(1,:) = conj(u_hat(end,:));

u = zeros(K,length(t));

for k = 1:K%K=5，要恢复5个分量
    u(k,:)=real(ifft(ifftshift(u_hat(:,k))));%复数的实部，快速傅里叶逆变换，逆零频平移（交换向量的左右两半部分）
end

%删除镜像部分 remove mirror part
u = u(:,T/4+1:3*T/4);%最终输出变量   分解模式的集合

% 再计算频谱recompute spectrum
clear u_hat;
for k = 1:K%K=5，要恢复5个分量
    u_hat(:,k)=fftshift(fft(u(k,:)))';%最终输出变量   模谱
end

end