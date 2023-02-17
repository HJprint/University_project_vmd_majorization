%% 基于遗传算法（Genetic Algorithm,GA)优化变分模态分解（variational mode decomposition,VMD）参数
clc
clear all
close all
tic%启动秒表计时器
% 读取数据
data=load('moni_noise.dat');%加载数据
s=data;%需要处理的信号
    
%% 设定遗传算法参数
maxgen=10;  % 进化代数，即迭代次数
sizepop=10; % 种群规模
pcross=0.8;  % 交叉概率选择(Pc:0-1)
pmutation=0.1; % 变异概率选择(Pm:0-1)
nvar=2; % 优化参数个数为2，分别为VMD的alpha和K
lenchrom=ones(1,nvar); % 个体长度 
% 参数范围（VMD的alpha和K）
bound=[100 5000;2 10];
% 种群初始化
individuals=struct('fitness',zeros(1,sizepop), 'chrom',[]);  %将种群信息定义为一个结构体
avgfitness=[]; %存储每一代种群的平均适应度
bestfitness=[]; %存储每一代种群的最佳适应度
bestchrom=[]; % 存储适应度最好的染色体
%初始化种群
for i=1:sizepop
    %随机产生一个种群
    individuals.chrom(i,:)=Code(lenchrom,bound);    %编码（binary和grey的编码结果为一个实数，float的编码结果为一个实数向量）
    x=individuals.chrom(i,:);
    %计算适应度
    individuals.fitness(i)=objfun(x,s);
end
%找最好的染色体
[bestfitness bestindex]=min(individuals.fitness);  %[m n]=min(b) m最小值 n列号
bestchrom=individuals.chrom(bestindex,:);  %最好的染色体
avgfitness=sum(individuals.fitness)/sizepop; %染色体的平均适应度
% 记录每一代进化中最好的适应度和平均适应度
trace=[avgfitness bestfitness]; 
%% 迭代寻优
start_time_train=cputime;
for i=1:maxgen
disp(['迭代次数：',num2str(i)])
% 选择
individuals=Select(individuals,sizepop); %选择后的新种群
avgfitness=sum(individuals.fitness)/sizepop;
% 交叉
individuals.chrom=Cross(pcross,lenchrom,individuals.chrom,sizepop,bound);
% 变异
individuals.chrom=Mutation(pmutation,lenchrom,individuals.chrom,sizepop,i,maxgen,bound);
% 计算适应度 
for j=1:sizepop
    x=individuals.chrom(j,:); %解码
    individuals.fitness(j)=objfun(x,s);
end
%找到最小和最大适应度的染色体及它们在种群中的位置
[newbestfitness,newbestindex]=min(individuals.fitness);
[worestfitness,worestindex]=max(individuals.fitness);
% 代替上一次进化中最好的染色体
if bestfitness>newbestfitness
    bestfitness=newbestfitness;
    bestchrom=individuals.chrom(newbestindex,:);
end
individuals.chrom(worestindex,:)=bestchrom;
individuals.fitness(worestindex)=bestfitness;
avgfitness=sum(individuals.fitness)/sizepop;
trace=[trace;avgfitness bestfitness]; %记录每一代进化中最好的适应度和平均适应度
end
%% 遗传算法结果分析 
figure
%解决中文字体显示问题
set(0,'defaultAxesFontName', 'Monospaced');
set(0,'defaultAxesFontSize', 10);

plot(trace(:,1),'r-o','linewidth',2)
hold on
xlabel('迭代次数')
ylabel('适应度')
legend('最佳适应度')
axis tight
grid on
%% 最优参数
xxx=bestchrom;
alpha=round(xxx(1));%需要优化的参数
K=round(xxx(2));%需要优化的参数
tau = 0;           
DC = 0;             
init = 1;          
tol = 1e-5;
%% VMD分解
[u, u_hat, omega] = VMD(s, alpha, tau, K, DC, init, tol);

figure;%建立幕布
for k=1:K%IMF
    subplot(K+3,1,k);plot(u(k,:),'k'); 
end
signal2=zeros(1,size(data,1) );%去噪后的信号
for k=1:K
    if  max(u(k,:))>10
        signal2=signal2+u(k,:);
    end
end
subplot(K+3,1,k+1);plot(data,'k');%原始信号 
subplot(K+3,1,k+2);plot(signal2,'k');


 toc%计算经过的时间


