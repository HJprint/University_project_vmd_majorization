%ͨ��WOA�㷨���ҵ����Ž��VMD
clear all;
clc;
tic%��������ʱ��
x=load('moni_noise.dat');%��������
%size(x)�ź�Ϊ1*2048
   signal=x;%����vmd�����õ��˲������˲���Ϊһ���źţ��ӱ��ж�ȡ�Ĳ���
   tau=0;%����vmd�����õ��˲���
   DC=0;%����vmd�����õ��˲���
   init=1;%����vmd�����õ��˲���
   tol=1e-5;%��1x10��-5�η�,����vmd�����õ��˲���
   
SearchAgents_no=10; % ��Ⱥ������Number of search agents
Max_iter=10; % ���������������ڼ���a,a2������zeros��������һ��1*10�ľ���Maximum numbef of iterations
dim=2; % ������Ҫ�Ż���������c��gλ������������zeros��������һ��1*2�ľ���number of your variables
lb=[2,100]; % ����ȡֵ�½磬��һ�����������������2��Ԫ��
ub=[10,5000]; % ����ȡֵ�Ͻ磬��һ�����������������2��Ԫ��
   
% The Whale Optimization Algorithm
% function [Leader_score,Leader_pos,Convergence_curve]=WOA(SearchAgents_no,Max_iter,lb,ub,dim,fobj)

% Ϊleader��ʼ��λ�������͵÷�,initialize position vector and score for the leader
Leader_pos=zeros(1,dim);%����һ��1*2�ľ���
Leader_score=inf; %infΪ�������+�ޣ�-infΪ����С��-��   change this to -inf for maximization problems


%��ʼ�����������λ��
Positions=initialization(SearchAgents_no,dim,ub,lb);%��ʼ��������PositionsΪһ��10*2��������ƾ���
%���initialization�������ص���PositionsΪһ��10*2��������ƾ���
%Positions��������10�����Ƶ�������Ϣ
%��������
Convergence_curve=zeros(1,Max_iter);%����һ��1*10��0����

t=0;% ѭ��������Loop counter

% ��ѭ��Main loop
while t<Max_iter%����������Ϊ10
    for i=1:size(Positions,1)%ѭ��ʮ��
        
        % ���س��������ռ�߽����������Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;%��������Ƚϣ���ÿһ����Ӧλ�õ�ֵ�Ƚ�
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        
        %����Ĵ���
        %����ÿ�����������Ŀ�꺯�� Calculate objective function for each search agent
        %round(Positions(i,2))�������ƣ�����ȡֵΪ�����㳤��1.5-2.0 ��
        %round(Positions(i,1))��Ҫ�ָ���ģʽ�������ֽ�ģ̬��IMF������
        %������������������
       [u, ~, ~] = VMD(signal,  round(Positions(i,2)), tau,  round(Positions(i,1)), DC, init, tol);
    
       %����ÿһ��position�ĵ÷�
     for ii=1:Positions(i,1)%size(Positions,1)=10��ÿһ��iʱ��ѭ��Positions��һ��������x��
        bao=hilbert(u(ii,:));
        bao=abs(bao);
        p=bao./sum(bao);
        e110(ii,:)=-sum(p.*log10(p));
     end
      
       fitness=min(e110);%����ÿһ��position�ĵ÷�
        
       % fitness=fobj(Positions(i,:));
        
        % �����쵼��Update the leader
        if fitness<Leader_score % Leader_score��ʼֵΪ�����Change this to > for maximization problem
            Leader_score=fitness; % Update alpha�����쵼�ߵĵ÷�
            Leader_pos=Positions(i,:);%10��λ��.�����쵼�ߵ�λ��
        end
        
    end
    
    % a decreases linearly fron 2 to 0 in Eq. (2.3)
    a=2-t*((2)/Max_iter); %�ڹ�ʽ(2.3)�У�a��2��0�����½���
    
    % a2 linearly dicreases from -1 to -2 to calculate t in Eq. (3.12)
    a2=-1+t*((-1)/Max_iter);%% a2��-1��-2���Եݼ������㹫ʽ(3.12)�е�t
    
    %���ݹ�ʽ���������������λ�� Update the Position of search agents 
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
                    rand_leader_index = floor(SearchAgents_no*rand()+1);%floor������ X ��ÿ��Ԫ���������뵽С�ڻ���ڸ�Ԫ�ص���ӽ�������
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
    
    %��������
    t=t+1;
    Convergence_curve(t)=Leader_score;%��������
     %[t Leader_score];ע�͵���
end     

%�õ������Ž�
bestc=Leader_pos(1,1);%����ķֽ�ģ̬��IMF������
bestg=Leader_pos(1,2);%����Ĵ�������
bestGWOaccuarcy=Leader_score;%���Ž�ĵ÷�

[u, ~, omega] = VMD(signal,  bestg, tau,  round(bestc), DC, init, tol);   

alpha=bestg;
K=round(bestc);%�� bestc ��ÿ��Ԫ����������Ϊ�����������

%�������������ʾ����
set(0,'defaultAxesFontName', 'Monospaced');
set(0,'defaultAxesFontSize', 10);
figure;%����Ļ��
for k=1:K%IMF
    subplot(K+3,1,k);plot(u(k,:),'k'); 
end
signal2=zeros(1,size(signal,1) );
for k=1:K
    if  max(u(k,:))>10
        signal2=signal2+u(k,:);
    end
end
subplot(K+3,1,k+1);plot(signal,'k');%ԭʼ�ź� 
subplot(K+3,1,k+2);plot(signal2,'k'); 

%��������
figure;
plot(Convergence_curve,'r-o','linewidth',2)
xlabel('��������')
ylabel('��Ӧ��')
legend('�����Ӧ��')



%out1=subplot(211);plot(x,'k'); xlabel('��������');ylabel('��ֵ');%xΪ�ź����ݣ�ԭʼ�ź�ͼ��k��ָ��ɫ��ɫ��x ��Ŀ̶ȷ�Χ�Ǵ� 1 �� length(x)��
%axis(out1,[0,length(signal),0,100]);%Ļ����x��y������

%out2=subplot(212);semilogy(Convergence_curve,'Color','k');xlabel('����������K��');ylabel('��Ӧ�Ⱥ���ֵ');%����ͼ
%axis(out2,[0,Max_iter,-inf,inf]);%Ļ����x��y������


toc%���㾭����ʱ��