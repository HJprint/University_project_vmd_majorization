% 天鹰优化器算法主程序
clear
close all
clc
tic%启动秒表计时器

%% I. 优化算法的参数设置
search_no=10; %种群大小
%适应度函数编号，设置范围: F1->F23
F_name='F1';
max_Iter=10;    %最大进化代数
%通过子函数读取适应度函数信息
signal=load('moni_noise.dat');%加载数据
[lb,ub,dim,objectf]=Get_FunctionDetails(F_name);
tau=0;%两处vmd函数用到此参数
DC=0;%两处vmd函数用到此参数
init=1;%两处vmd函数用到此参数
tol=1e-5;%是1x10的-5次方,两处vmd函数用到此参数

%% II. 通过优化算法对适应度函数寻优
[best_position,best_score,curve]=AO(signal,tau,DC,init,tol,search_no,max_Iter,lb,ub,dim,objectf); 

%% III. 绘制进化曲线
alpha=best_position(1,2);
K=round(best_position(1,1));
[u, ~, omega] = VMD(signal,  alpha, tau,  K, DC, init, tol); 
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
subplot(K+3,1,k+1);plot(signal,'k');%原始信号 
subplot(K+3,1,k+2);plot(signal2,'k');


%收敛曲线
figure;
plot(curve,'r-o','linewidth',2)
xlabel('迭代次数')
ylabel('适应度')
legend('最佳适应度')

toc%计算经过的时间




