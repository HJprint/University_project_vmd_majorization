% 主程序 GWO灰狼优化算法
clear
close all
clc
tic%启动秒表计时器
%vmd参数初始化
 x=load('moni_noise.dat');%加载数据
 signal=x;%两处vmd函数用到此参数，此参数为一个信号，从表中读取的参数
   tau=0;%两处vmd函数用到此参数
   DC=0;%两处vmd函数用到此参数
   init=1;%两处vmd函数用到此参数
   tol=1e-5;%是1x10的-5次方,两处vmd函数用到此参数
 
   %GWO参数初始化
SearchAgents_no = 10 ; %　种群规模
dim = 2 ; % 粒子维度.两个优化参数
Max_iter = 10; %　迭代次数
lb=[2,100]; % 参数取值下界
ub=[10,5000]; % 参数取值上界
 
% 初始化三匹头狼的位置
Alpha_pos=zeros(1,dim);%第一头狼的位置
Alpha_score=inf; %第一头狼的得分
 
Beta_pos=zeros(1,dim);
Beta_score=inf; 
 
Delta_pos=zeros(1,dim);
Delta_score=inf; 
 
%初始化种群的位置
Positions = initialization(SearchAgents_no,dim,ub,lb); % 初始化粒子群的位置
 %收敛曲线
Convergence_curve = zeros(Max_iter,1);
 
% 开始循环
for l=1:Max_iter
    for i=1:size(Positions,1)  
        
       % 返回超出搜索空间边界的搜索代理Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;               
        
      %计算每一个position的得分
        %fitness=sum(Positions(i,:).^2);
        [u, ~, ~] = VMD(signal,  round(Positions(i,2)), tau,  round(Positions(i,1)), DC, init, tol);
     for ii=1:Positions(i,1)%size(Positions,1)=10，每一个i时，循环Positions第一个参数的x次
        bao=hilbert(u(ii,:));
        bao=abs(bao);
        p=bao./sum(bao);
        e110(ii,:)=-sum(p.*log10(p));
     end
       fitness=min(e110);%计算每一个position的得分
       
        % Update Alpha, Beta, and Delta
        if fitness<Alpha_score %小于A狼------A
            Alpha_score=fitness; % Update alpha
            Alpha_pos=Positions(i,:);
        end
        
        if fitness>Alpha_score && fitness<Beta_score%A，B狼之间 -------B
            Beta_score=fitness; % Update beta
            Beta_pos=Positions(i,:);
        end
        
        if fitness>Alpha_score && fitness>Beta_score && fitness<Delta_score %B,D之间-------D
            Delta_score=fitness; % Update delta
            Delta_pos=Positions(i,:);
        end
    end
    
    
    a=2-l*((2)/Max_iter); % A从2线性减少到0  a decreases linearly fron 2 to 0
    
    % Update the Position of search agents including omegas
    for i=1:size(Positions,1)
        for j=1:size(Positions,2)     
                       
            r1=rand(); % r1 is a random number in [0,1]
            r2=rand(); % r2 is a random number in [0,1]
            
            A1=2*a*r1-a; % Equation (3.3)
            C1=2*r2; % Equation (3.4)
            
            D_alpha=abs(C1*Alpha_pos(j)-Positions(i,j)); % Equation (3.5)-part 1
            X1=Alpha_pos(j)-A1*D_alpha; % Equation (3.6)-part 1
                       
            r1=rand();
            r2=rand();
            
            A2=2*a*r1-a; % Equation (3.3)
            C2=2*r2; % Equation (3.4)
            
            D_beta=abs(C2*Beta_pos(j)-Positions(i,j)); % Equation (3.5)-part 2
            X2=Beta_pos(j)-A2*D_beta; % Equation (3.6)-part 2       
            
            r1=rand();
            r2=rand(); 
            
            A3=2*a*r1-a; % Equation (3.3)
            C3=2*r2; % Equation (3.4)
            
            D_delta=abs(C3*Delta_pos(j)-Positions(i,j)); % Equation (3.5)-part 3
            X3=Delta_pos(j)-A3*D_delta; % Equation (3.5)-part 3             
            
            Positions(i,j)=(X1+X2+X3)/3;% Equation (3.7)
            
        end
    end
  
    Convergence_curve(l)=Alpha_score;%收敛曲线
    %disp(['Iteration = ' num2str(l)  ', Evaluations = ' num2str(Alpha_score)]);
 
end
%得到的最优解
bestc=Alpha_pos(1,1);%输入的分解模态（IMF）个数
bestg=Alpha_pos(1,2);%输入的带宽限制
bestGWOaccuarcy=Alpha_score;%最优解的得分
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
signal2=zeros(1,size(signal,1) );%去噪后的信号
for k=1:K
    if  max(u(k,:))>10
        signal2=signal2+u(k,:);
    end
end
subplot(K+3,1,k+1);plot(x,'k');%原始信号 
subplot(K+3,1,k+2);plot(signal2,'k');

%收敛曲线
figure;
plot(Convergence_curve,'r-o','linewidth',2)
xlabel('迭代次数')
ylabel('适应度')
legend('最佳适应度')

 toc%计算经过的时间