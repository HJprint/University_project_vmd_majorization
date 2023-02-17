% ��ӥ�Ż����㷨������
clear
close all
clc
tic%��������ʱ��

%% I. �Ż��㷨�Ĳ�������
search_no=10; %��Ⱥ��С
%��Ӧ�Ⱥ�����ţ����÷�Χ: F1->F23
F_name='F1';
max_Iter=10;    %����������
%ͨ���Ӻ�����ȡ��Ӧ�Ⱥ�����Ϣ
signal=load('moni_noise.dat');%��������
[lb,ub,dim,objectf]=Get_FunctionDetails(F_name);
tau=0;%����vmd�����õ��˲���
DC=0;%����vmd�����õ��˲���
init=1;%����vmd�����õ��˲���
tol=1e-5;%��1x10��-5�η�,����vmd�����õ��˲���

%% II. ͨ���Ż��㷨����Ӧ�Ⱥ���Ѱ��
[best_position,best_score,curve]=AO(signal,tau,DC,init,tol,search_no,max_Iter,lb,ub,dim,objectf); 

%% III. ���ƽ�������
alpha=best_position(1,2);
K=round(best_position(1,1));
[u, ~, omega] = VMD(signal,  alpha, tau,  K, DC, init, tol); 
figure;%����Ļ��
for k=1:K%IMF
    subplot(K+3,1,k);plot(u(k,:),'k'); 
end
signal2=zeros(1,size(signal,1) );%ȥ�����ź�
for k=1:K
    if  max(u(k,:))>10
        signal2=signal2+u(k,:);
    end
end
subplot(K+3,1,k+1);plot(signal,'k');%ԭʼ�ź� 
subplot(K+3,1,k+2);plot(signal2,'k');


%��������
figure;
plot(curve,'r-o','linewidth',2)
xlabel('��������')
ylabel('��Ӧ��')
legend('�����Ӧ��')

toc%���㾭����ʱ��




