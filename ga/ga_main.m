%% �����Ŵ��㷨��Genetic Algorithm,GA)�Ż����ģ̬�ֽ⣨variational mode decomposition,VMD������
clc
clear all
close all
tic%��������ʱ��
% ��ȡ����
data=load('moni_noise.dat');%��������
s=data;%��Ҫ������ź�
    
%% �趨�Ŵ��㷨����
maxgen=10;  % ��������������������
sizepop=10; % ��Ⱥ��ģ
pcross=0.8;  % �������ѡ��(Pc:0-1)
pmutation=0.1; % �������ѡ��(Pm:0-1)
nvar=2; % �Ż���������Ϊ2���ֱ�ΪVMD��alpha��K
lenchrom=ones(1,nvar); % ���峤�� 
% ������Χ��VMD��alpha��K��
bound=[100 5000;2 10];
% ��Ⱥ��ʼ��
individuals=struct('fitness',zeros(1,sizepop), 'chrom',[]);  %����Ⱥ��Ϣ����Ϊһ���ṹ��
avgfitness=[]; %�洢ÿһ����Ⱥ��ƽ����Ӧ��
bestfitness=[]; %�洢ÿһ����Ⱥ�������Ӧ��
bestchrom=[]; % �洢��Ӧ����õ�Ⱦɫ��
%��ʼ����Ⱥ
for i=1:sizepop
    %�������һ����Ⱥ
    individuals.chrom(i,:)=Code(lenchrom,bound);    %���루binary��grey�ı�����Ϊһ��ʵ����float�ı�����Ϊһ��ʵ��������
    x=individuals.chrom(i,:);
    %������Ӧ��
    individuals.fitness(i)=objfun(x,s);
end
%����õ�Ⱦɫ��
[bestfitness bestindex]=min(individuals.fitness);  %[m n]=min(b) m��Сֵ n�к�
bestchrom=individuals.chrom(bestindex,:);  %��õ�Ⱦɫ��
avgfitness=sum(individuals.fitness)/sizepop; %Ⱦɫ���ƽ����Ӧ��
% ��¼ÿһ����������õ���Ӧ�Ⱥ�ƽ����Ӧ��
trace=[avgfitness bestfitness]; 
%% ����Ѱ��
start_time_train=cputime;
for i=1:maxgen
disp(['����������',num2str(i)])
% ѡ��
individuals=Select(individuals,sizepop); %ѡ��������Ⱥ
avgfitness=sum(individuals.fitness)/sizepop;
% ����
individuals.chrom=Cross(pcross,lenchrom,individuals.chrom,sizepop,bound);
% ����
individuals.chrom=Mutation(pmutation,lenchrom,individuals.chrom,sizepop,i,maxgen,bound);
% ������Ӧ�� 
for j=1:sizepop
    x=individuals.chrom(j,:); %����
    individuals.fitness(j)=objfun(x,s);
end
%�ҵ���С�������Ӧ�ȵ�Ⱦɫ�弰��������Ⱥ�е�λ��
[newbestfitness,newbestindex]=min(individuals.fitness);
[worestfitness,worestindex]=max(individuals.fitness);
% ������һ�ν�������õ�Ⱦɫ��
if bestfitness>newbestfitness
    bestfitness=newbestfitness;
    bestchrom=individuals.chrom(newbestindex,:);
end
individuals.chrom(worestindex,:)=bestchrom;
individuals.fitness(worestindex)=bestfitness;
avgfitness=sum(individuals.fitness)/sizepop;
trace=[trace;avgfitness bestfitness]; %��¼ÿһ����������õ���Ӧ�Ⱥ�ƽ����Ӧ��
end
%% �Ŵ��㷨������� 
figure
%�������������ʾ����
set(0,'defaultAxesFontName', 'Monospaced');
set(0,'defaultAxesFontSize', 10);

plot(trace(:,1),'r-o','linewidth',2)
hold on
xlabel('��������')
ylabel('��Ӧ��')
legend('�����Ӧ��')
axis tight
grid on
%% ���Ų���
xxx=bestchrom;
alpha=round(xxx(1));%��Ҫ�Ż��Ĳ���
K=round(xxx(2));%��Ҫ�Ż��Ĳ���
tau = 0;           
DC = 0;             
init = 1;          
tol = 1e-5;
%% VMD�ֽ�
[u, u_hat, omega] = VMD(s, alpha, tau, K, DC, init, tol);

figure;%����Ļ��
for k=1:K%IMF
    subplot(K+3,1,k);plot(u(k,:),'k'); 
end
signal2=zeros(1,size(data,1) );%ȥ�����ź�
for k=1:K
    if  max(u(k,:))>10
        signal2=signal2+u(k,:);
    end
end
subplot(K+3,1,k+1);plot(data,'k');%ԭʼ�ź� 
subplot(K+3,1,k+2);plot(signal2,'k');


 toc%���㾭����ʱ��


