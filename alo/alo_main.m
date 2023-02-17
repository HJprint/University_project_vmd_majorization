%___________________________________________________________________%
%  Ant Lion Optimizer (ALO) source codes demo version 1.0           %
%    蚁狮优化算法                                                               %
%  Developed in MATLAB R2011b(7.13)                                 %
%                                                                   %
%  Author and programmer: Seyedali Mirjalili                        %
%                                                                   %
%         e-Mail: ali.mirjalili@gmail.com                           %
%                 seyedali.mirjalili@griffithuni.edu.au             %
%                                                                   %
%       Homepage: http://www.alimirjalili.com                       %
%                                                                   %
%   Main paper:                                                     %
%                                                                   %
%   S. Mirjalili, The Ant Lion Optimizer                            %
%   Advances in Engineering Software , in press,2015                %
%   DOI: http://dx.doi.org/10.1016/j.advengsoft.2015.01.010         %
%                                                                   %
%___________________________________________________________________%

% You can simply define your cost in a seperate file and load its handle to fobj 
% The initial parameters that you need are:
%__________________________________________
% fobj = @YourCostFunction
% dim = number of your variables
% Max_iteration = maximum number of generations
% SearchAgents_no = number of search agents
% lb=[lb1,lb2,...,lbn] where lbn is the lower bound of variable n
% ub=[ub1,ub2,...,ubn] where ubn is the upper bound of variable n
% If all the variables have equal lower bound you can just
% define lb and ub as two single number numbers

% To run ALO: [Best_score,Best_pos,cg_curve]=ALO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj)
%__________________________________________

clear all 
clc
tic%启动秒表计时器
x=load('moni_noise.dat');%加载数据
%size(x)信号为1*2048
signal=x;%两处vmd函数用到此参数，此参数为一个信号，从表中读取的参数
tau=0;%两处vmd函数用到此参数
DC=0;%两处vmd函数用到此参数
init=1;%两处vmd函数用到此参数
tol=1e-5;%是1x10的-5次方,两处vmd函数用到此参数
   
SearchAgents_no=10; % Number of search agents
Function_name='F1'; % Name of the test function that can be from F1 to F23 (Table 1,2,3 in the paper)
Max_iteration=10; % Maximum numbef of iterations
% Load details of the selected benchmark function
[lb,ub,dim,fobj]=Get_Functions_details(Function_name);

[Best_score,Best_pos,cg_curve]=ALO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj,signal,tau,DC,init,tol);

bestc=Best_pos(1,1);%输入的分解模态（IMF）个数
bestg=Best_pos(1,2);%输入的带宽限制
bestGWOaccuarcy=Best_score;%最优解的得分

[u, ~, omega] = VMD(signal,  bestg, tau,  round(bestc), DC, init, tol);   

alpha=bestg;
K=round(bestc);%将 bestc 的每个元素四舍五入为最近的整数。

figure;%建立幕布
for k=1:K%IMF
    subplot(K+3,1,k);plot(u(k,:),'k'); 
end
signal2=zeros(1,size(x,1) );%去噪后的信号
for k=1:K
    if  max(u(k,:))>10
        signal2=signal2+u(k,:);
    end
end
subplot(K+3,1,k+1);plot(x,'k');%原始信号 
subplot(K+3,1,k+2);plot(signal2,'k');%去噪后的信号

%收敛曲线
figure;
plot(cg_curve,'r-o','linewidth',2)
xlabel('迭代次数')
ylabel('适应度')
legend('最佳适应度')
        
toc%计算经过的时间


