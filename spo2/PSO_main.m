%% ��ջ���
clear;
clc;
tic%��������ʱ��

%% ��������
%signal=load('moni_noise.dat');%��������
[signal1,fs]=audioread('he.mp3');%��������
signal=signal1(:,1);
figure
plot(signal1)
w=0.9;%Ȩֵ ��Ӱ��PSO ��ȫ����ֲ����������� ֵ�ϴ�ȫ����������ǿ���ֲ�����������;��֮����ֲ�����������ǿ����ȫ����������������
c1=0.1;%���ٶȣ�Ӱ�������ٶ�
c2=0.1;
dim=2;%6ά
swarmsize=10;%����Ⱥ��ģ����ʾ��10����Ŀռ�
maxiter=10;%���ѭ��������Ӱ��ʱ��
minfit=0.001;%��С��Ӧֵ
vmax=0.01;
vmin=-0.01;
ub=[10,5000];%���������������
lb=[2,100];%����������С����
tau=0;%����vm d�����õ��˲���
DC=0;%����vmd�����õ��˲���
init=1;%����vmd�����õ��˲���
tol=1e-5;%��1x10��-5�η�,����vmd�����õ��˲���

%% ��Ⱥ��ʼ��
swarm=initialization(swarmsize,dim,ub,lb);
vstep=rand(swarmsize,dim)*(vmax-vmin)+vmin;%����Ⱥ�ٶȾ���
fswarm=zeros(swarmsize,1);%Ԥ��վ��󣬴����Ӧֵ


%% ���弫ֵ��Ⱥ�弫ֵ
[bestf,bestindex]=min(fswarm);%�����Ӧֵ�е���С��Ӧֵ���ͣ������ڵ�����
gbest=swarm;%��ʱ�ĸ������Ž�Ϊ�Լ�
fgbest=fswarm;%��ʱ�ĸ���������Ӧֵ
zbest=swarm(bestindex,:);%�������еĶ�Ӧ�Ľ�������У�ȫ����ѽ�
fzbest=bestf;%ȫ��������Ӧֵ


%% ����Ѱ��
iter=0;
yfitness=zeros(1,maxiter);%1��10�о��󣬴��10������ֵ�Ŀռ����

while((iter<maxiter)&&(fzbest>minfit))%����10��
    for j=1:swarmsize%ѭ��10������
        % �ٶȸ���
        vstep(j,:)=w*vstep(j,:)+c1*rand*(gbest(j,:)-swarm(j,:))+c2*rand*(zbest-swarm(j,:));
        if vstep(j,:)>vmax  
            vstep(j,:)=vmax;%�ٶ�����
        end
        if vstep(j,:)<vmin
            vstep(j,:)=vmin;
        end
        % λ�ø���
        swarm(j,:)=swarm(j,:)+vstep(j,:);
        for k=1:dim
            if swarm(j,k)>ub(k)
                swarm(j,k)=ub(k);%λ������
            end
            if swarm(j,k)<lb(k)
                swarm(j,k)=lb(k);
            end
        end
       
        
        % ���س��������ռ�߽����������
        Flag4ub=swarm(j,:)>ub;%��������Ƚϣ���ÿһ����Ӧλ�õ�ֵ�Ƚ�
        Flag4lb=swarm(j,:)<lb;
        swarm(j,:)=(swarm(j,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
         % ��Ӧֵ
         for ii=1:swarm(j,1)%size(Positions,1)=10��ÿһ��iʱ��ѭ��Positions��һ��������x��
            bao=hilbert(u(ii,:));
            bao=abs(bao);
            p=bao./sum(bao);
            e110(ii)=-sum(p.*log10(p))
         end
       fswarm(j,:)=min(e110);%����ÿһ��position�ĵ÷�
       %x=swarm(j,:); 
       %fswarm(j,:)=objfun(x,signal);%����ÿһ��position�ĵ÷�
        % ���ڴ˴�����Լ��������������Լ���������������Ӧֵ����
        
        %
        % �������Ÿ���
        if fswarm(j)<fgbest(j) %�����ǰ�ĺ���ֵ�ȸ�������ֵС
            gbest(j,:)=swarm(j,:);%�������Ž����
            fgbest(j)=fswarm(j);%��������ֵ����
        end
        % Ⱥ�����Ÿ���
        if fswarm(j)<fzbest%�����ǰ�ĺ���ֵ��Ⱥ������ֵ��
            zbest=swarm(j,:);%Ⱥ�����Ž����
            fzbest=fswarm(j);%Ⱥ������ֵ����
        end
    end
    iter=iter+1;
    yfitness(iter)=fzbest;

end


%% ��ͼ
%% ���Ų���
xxx=zbest;
alpha=round(xxx(2));%��Ҫ�Ż��Ĳ���
K=round(xxx(1));%��Ҫ�Ż��Ĳ���
%% VMD�ֽ�
[u, u_hat, omega] = VMD(signal, alpha, tau, K, DC, init, tol);
omega;

figure(1);%����Ļ��
%�������������ʾ����
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
subplot(2,1,1);plot(signal,'k');%ԭʼ�ź� 
subplot(2,1,2);plot(signal2,'k');%ȥ���ź� 

figure(3);
plot(yfitness,'r-o','linewidth',2)
xlabel('��������')
ylabel('��Ӧ��')
legend('�����Ӧ��')

%axis tight%�Զ�����x���y��ķ�Χʹͼ����������ռ��������ʾ�ռ�


 toc%���㾭����ʱ��







