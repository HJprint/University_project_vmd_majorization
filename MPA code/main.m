%_________________________________________________________________________
%  Marine Predators Algorithm source code (Developed in MATLAB R2015a)
%
%  programming: Afshin Faramarzi & Seyedali Mirjalili
%
% paper:
%  A. Faramarzi, M. Heidarinejad, S. Mirjalili, A.H. Gandomi, 
%  Marine Predators Algorithm: A Nature-inspired Metaheuristic
%  Expert Systems with Applications
%  DOI: doi.org/10.1016/j.eswa.2020.113377
%  
%  E-mails: afaramar@hawk.iit.edu            (Afshin Faramarzi)
%           muh182@iit.edu                   (Mohammad Heidarinejad)
%           ali.mirjalili@laureate.edu.au    (Seyedali Mirjalili) 
%           gandomi@uts.edu.au               (Amir H Gandomi)
%_________________________________________________________________________

% --------------------------------------------
% fobj = @YourCostFunction
% dim = number of your variables
% Max_iteration = maximum number of iterations
% SearchAgents_no = number of search agents
% lb=[lb1,lb2,...,lbn] where lbn is the lower bound of variable n
% ub=[ub1,ub2,...,ubn] where ubn is the upper bound of variable n
% ---------------------------------------------------------

clear all
clc
tic%启动秒表计时器
format long

x=load('moni_noise.dat');%加载数据
%size(x)信号为1*2048
   signal=x;%两处vmd函数用到此参数，此参数为一个信号，从表中读取的参数
   tau=0;%两处vmd函数用到此参数
   DC=0;%两处vmd函数用到此参数
   init=1;%两处vmd函数用到此参数
   tol=1e-5;%是1x10的-5次方,两处vmd函数用到此参数

SearchAgents_no=10; % 种群规模Number of search agents

Function_name='F1';
   
Max_iteration=10; % 迭代次数Maximum number of iterations

[lb,ub,dim,fobj]=Get_Functions_details(Function_name);

[Best_score,Best_pos,Convergence_curve]=MPA(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);

%得到的最优解
bestc=Best_pos(1,1);%输入的分解模态（IMF）个数
bestg=Best_pos(1,2);%输入的带宽限制
bestGWOaccuarcy=Best_score;%最优解的得分

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


 toc%计算经过的时间
