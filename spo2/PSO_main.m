%% 清空环境
clear;
clc;
tic%启动秒表计时器

%% 参数设置
%signal=load('moni_noise.dat');%加载数据
[signal1,fs]=audioread('he.mp3');%加载数据
signal=signal1(:,1);
figure
plot(signal1)
w=0.9;%权值 将影响PSO 的全局与局部搜优能力， 值较大，全局搜优能力强，局部搜优能力弱;反之，则局部搜优能力增强，而全局搜优能力减弱。
c1=0.1;%加速度，影响收敛速度
c2=0.1;
dim=2;%6维
swarmsize=10;%粒子群规模，表示有10个解的空间
maxiter=10;%最大循环次数，影响时间
minfit=0.001;%最小适应值
vmax=0.01;
vmin=-0.01;
ub=[10,5000];%解向量的最大限制
lb=[2,100];%解向量的最小限制
tau=0;%两处vm d函数用到此参数
DC=0;%两处vmd函数用到此参数
init=1;%两处vmd函数用到此参数
tol=1e-5;%是1x10的-5次方,两处vmd函数用到此参数

%% 种群初始化
swarm=initialization(swarmsize,dim,ub,lb);
vstep=rand(swarmsize,dim)*(vmax-vmin)+vmin;%粒子群速度矩阵
fswarm=zeros(swarmsize,1);%预设空矩阵，存放适应值


%% 个体极值和群体极值
[bestf,bestindex]=min(fswarm);%求得适应值中的最小适应值，和，其所在的序列
gbest=swarm;%暂时的个体最优解为自己
fgbest=fswarm;%暂时的个体最优适应值
zbest=swarm(bestindex,:);%所在序列的对应的解矩阵序列，全局最佳解
fzbest=bestf;%全局最优适应值


%% 迭代寻优
iter=0;
yfitness=zeros(1,maxiter);%1行10列矩阵，存放10个最优值的空间矩阵

while((iter<maxiter)&&(fzbest>minfit))%迭代10次
    for j=1:swarmsize%循环10个粒子
        % 速度更新
        vstep(j,:)=w*vstep(j,:)+c1*rand*(gbest(j,:)-swarm(j,:))+c2*rand*(zbest-swarm(j,:));
        if vstep(j,:)>vmax  
            vstep(j,:)=vmax;%速度限制
        end
        if vstep(j,:)<vmin
            vstep(j,:)=vmin;
        end
        % 位置更新
        swarm(j,:)=swarm(j,:)+vstep(j,:);
        for k=1:dim
            if swarm(j,k)>ub(k)
                swarm(j,k)=ub(k);%位置限制
            end
            if swarm(j,k)<lb(k)
                swarm(j,k)=lb(k);
            end
        end
       
        
        % 返回超出搜索空间边界的搜索代理
        Flag4ub=swarm(j,:)>ub;%两个矩阵比较，是每一个对应位置的值比较
        Flag4lb=swarm(j,:)<lb;
        swarm(j,:)=(swarm(j,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
         % 适应值
         for ii=1:swarm(j,1)%size(Positions,1)=10，每一个i时，循环Positions第一个参数的x次
            bao=hilbert(u(ii,:));
            bao=abs(bao);
            p=bao./sum(bao);
            e110(ii)=-sum(p.*log10(p))
         end
       fswarm(j,:)=min(e110);%计算每一个position的得分
       %x=swarm(j,:); 
       %fswarm(j,:)=objfun(x,signal);%计算每一个position的得分
        % 可在此处增加约束条件，若满足约束条件，则进行适应值计算
        
        %
        % 个体最优更新
        if fswarm(j)<fgbest(j) %如果当前的函数值比个体最优值小
            gbest(j,:)=swarm(j,:);%个体最优解更新
            fgbest(j)=fswarm(j);%个体最优值更新
        end
        % 群体最优更新
        if fswarm(j)<fzbest%如果当前的函数值比群体最优值大
            zbest=swarm(j,:);%群体最优解更新
            fzbest=fswarm(j);%群体最优值更新
        end
    end
    iter=iter+1;
    yfitness(iter)=fzbest;

end


%% 画图
%% 最优参数
xxx=zbest;
alpha=round(xxx(2));%需要优化的参数
K=round(xxx(1));%需要优化的参数
%% VMD分解
[u, u_hat, omega] = VMD(signal, alpha, tau, K, DC, init, tol);
omega;

figure(1);%建立幕布
%解决中文字体显示问题
set(0,'defaultAxesFontName', 'Monospaced');
set(0,'defaultAxesFontSize', 10);
for k=1:K%IMF
    subplot(K,1,k);plot(u(k,:),'k'); 
end
signal2=zeros(1,size(signal,1) );
for k=1:K
    if  max(u(k,:))>0.2
        signal2=signal2+u(k,:);
    end
end

fname='ok.wav';
audiowrite(fname,signal2,fs);
[signal2,fs]=audioread(fname);
sound(signal2,fs);

figure(2);
subplot(2,1,1);plot(signal,'k');%原始信号 
subplot(2,1,2);plot(signal2,'k');%去噪信号 

figure(3);
plot(yfitness,'r-o','linewidth',2)
xlabel('迭代次数')
ylabel('适应度')
legend('最佳适应度')

%axis tight%自动设置x轴和y轴的范围使图形区域正好占满整个显示空间


 toc%计算经过的时间







