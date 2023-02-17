%通过WOA算法，找到最优解给VMD
clear all;
clc;
tic%启动秒表计时器
x=load('moni_noise.dat');%加载数据
%size(x)信号为1*2048
   signal=x;%两处vmd函数用到此参数，此参数为一个信号，从表中读取的参数
   tau=0;%两处vmd函数用到此参数
   DC=0;%两处vmd函数用到此参数
   init=1;%两处vmd函数用到此参数
   tol=1e-5;%是1x10的-5次方,两处vmd函数用到此参数
   
SearchAgents_no=10; % 种群数量，Number of search agents
Max_iter=10; % 最大迭代次数，用于计算a,a2。用于zeros函数生成一个1*10的矩阵。Maximum numbef of iterations
dim=2; % 此例需要优化两个参数c和g位置向量，用于zeros函数生成一个1*2的矩阵。number of your variables
lb=[2,100]; % 参数取值下界，是一个向量，这个向量有2个元素
ub=[10,5000]; % 参数取值上界，是一个向量，这个向量有2个元素
   
% The Whale Optimization Algorithm
% function [Leader_score,Leader_pos,Convergence_curve]=WOA(SearchAgents_no,Max_iter,lb,ub,dim,fobj)

% 为leader初始化位置向量和得分,initialize position vector and score for the leader
Leader_pos=zeros(1,dim);%生成一个1*2的矩阵。
Leader_score=inf; %inf为无穷大量+∞，-inf为无穷小量-∞   change this to -inf for maximization problems


%初始化搜索代理的位置
Positions=initialization(SearchAgents_no,dim,ub,lb);%初始化函数，Positions为一个10*2的随机限制矩阵
%最后initialization函数返回的是Positions为一个10*2的随机限制矩阵
%Positions里面存放着10个限制的坐标信息
%收敛曲线
Convergence_curve=zeros(1,Max_iter);%生成一个1*10的0矩阵。

t=0;% 循环计数器Loop counter

% 主循环Main loop
while t<Max_iter%最大迭代次数为10
    for i=1:size(Positions,1)%循环十次
        
        % 返回超出搜索空间边界的搜索代理Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;%两个矩阵比较，是每一个对应位置的值比较
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        
        %多出的代码
        %计算每个搜索代理的目标函数 Calculate objective function for each search agent
        %round(Positions(i,2))带宽限制，经验取值为抽样点长度1.5-2.0 倍
        %round(Positions(i,1))需要恢复的模式个数，分解模态（IMF）个数
        %不接收那两个参数了
       [u, ~, ~] = VMD(signal,  round(Positions(i,2)), tau,  round(Positions(i,1)), DC, init, tol);
    
       %计算每一个position的得分
     for ii=1:Positions(i,1)%size(Positions,1)=10，每一个i时，循环Positions第一个参数的x次
        bao=hilbert(u(ii,:));
        bao=abs(bao);
        p=bao./sum(bao);
        e110(ii,:)=-sum(p.*log10(p));
     end
      
       fitness=min(e110);%计算每一个position的得分
        
       % fitness=fobj(Positions(i,:));
        
        % 更新领导者Update the leader
        if fitness<Leader_score % Leader_score初始值为无穷大Change this to > for maximization problem
            Leader_score=fitness; % Update alpha更新领导者的得分
            Leader_pos=Positions(i,:);%10个位置.更新领导者的位置
        end
        
    end
    
    % a decreases linearly fron 2 to 0 in Eq. (2.3)
    a=2-t*((2)/Max_iter); %在公式(2.3)中，a从2到0线性下降。
    
    % a2 linearly dicreases from -1 to -2 to calculate t in Eq. (3.12)
    a2=-1+t*((-1)/Max_iter);%% a2从-1到-2线性递减来计算公式(3.12)中的t
    
    %根据公式，更新搜索代理的位置 Update the Position of search agents 
    for i=1:size(Positions,1)%size(Positions,1)=10
        r1=rand(); % r1 is a random number in [0,1]
        r2=rand(); % r2 is a random number in [0,1]
        
        A=2*a*r1-a;  % Eq. (2.3) in the paper
        C=2*r2;      % Eq. (2.4) in the paper
        
        
        b=1;               %  parameters in Eq. (2.5)
        l=(a2-1)*rand+1;   %  parameters in Eq. (2.5)
        
        p = rand();        % p in Eq. (2.6)
        
        for j=1:size(Positions,2)%size(Positions,2)=2
            
            if p<0.5   
                if abs(A)>=1
                    rand_leader_index = floor(SearchAgents_no*rand()+1);%floor函数将 X 的每个元素四舍五入到小于或等于该元素的最接近整数。
                    X_rand = Positions(rand_leader_index, :);
                    D_X_rand=abs(C*X_rand(j)-Positions(i,j)); % Eq. (2.7)
                    Positions(i,j)=X_rand(j)-A*D_X_rand;      % Eq. (2.8)
                    
                elseif abs(A)<1
                    D_Leader=abs(C*Leader_pos(j)-Positions(i,j)); % Eq. (2.1)
                    Positions(i,j)=Leader_pos(j)-A*D_Leader;      % Eq. (2.2)
                end
                
            elseif p>=0.5
              
                distance2Leader=abs(Leader_pos(j)-Positions(i,j));
                % Eq. (2.5)
                Positions(i,j)=distance2Leader*exp(b.*l).*cos(l.*2*pi)+Leader_pos(j);
                
            end
            
        end
    end
    
    %收敛曲线
    t=t+1;
    Convergence_curve(t)=Leader_score;%收敛曲线
     %[t Leader_score];注释掉了
end     

%得到的最优解
bestc=Leader_pos(1,1);%输入的分解模态（IMF）个数
bestg=Leader_pos(1,2);%输入的带宽限制
bestGWOaccuarcy=Leader_score;%最优解的得分

[u, ~, omega] = VMD(signal,  bestg, tau,  round(bestc), DC, init, tol);   

alpha=bestg;
K=round(bestc);%将 bestc 的每个元素四舍五入为最近的整数。

%解决中文字体显示问题
set(0,'defaultAxesFontName', 'Monospaced');
set(0,'defaultAxesFontSize', 10);
figure;%建立幕布
for k=1:K%IMF
    subplot(K+3,1,k);plot(u(k,:),'k'); 
end
signal2=zeros(1,size(signal,1) );
for k=1:K
    if  max(u(k,:))>10
        signal2=signal2+u(k,:);
    end
end
subplot(K+3,1,k+1);plot(signal,'k');%原始信号 
subplot(K+3,1,k+2);plot(signal2,'k'); 

%收敛曲线
figure;
plot(Convergence_curve,'r-o','linewidth',2)
xlabel('迭代次数')
ylabel('适应度')
legend('最佳适应度')



%out1=subplot(211);plot(x,'k'); xlabel('采样点数');ylabel('幅值');%x为信号数据，原始信号图，k是指颜色黑色，x 轴的刻度范围是从 1 至 length(x)。
%axis(out1,[0,length(signal),0,100]);%幕布的x，y轴限制

%out2=subplot(212);semilogy(Convergence_curve,'Color','k');xlabel('迭代次数（K）');ylabel('适应度函数值');%迭代图
%axis(out2,[0,Max_iter,-inf,inf]);%幕布的x，y轴限制


toc%计算经过的时间