%_________________________________________________________________________%
%  Whale Optimization Algorithm (WOA) source codes demo 1.0               %
%                                                                         %
%  Developed in MATLAB R2011b(7.13)                                       %
%                                                                         %
%  Author and programmer: Seyedali Mirjalili                              %
%                                                                         %
%         e-Mail: ali.mirjalili@gmail.com                                 %
%                 seyedali.mirjalili@griffithuni.edu.au                   %
%                                                                         %
%       Homepage: http://www.alimirjalili.com                             %
%                                                                         %
%   Main paper: S. Mirjalili, A. Lewis                                    %
%               The Whale Optimization Algorithm,                         %
%               Advances in Engineering Software , in press,              %
%               DOI: http://dx.doi.org/10.1016/j.advengsoft.2016.01.008   %
%                                                                         %
%_________________________________________________________________________%

function [out1, out2, out3] = WOAVMD(x) 
   signal=x;
   tau=0;
   DC=0;
   init=1;
   tol=1e-5;
   
SearchAgents_no=10; % 种群数量，Number of search agents
Max_iter=10; % 最大迭代次数，Maximum numbef of iterations
dim=2; % 此例需要优化两个参数c和g，number of your variables
lb=[2,200]; % 参数取值下界
ub=[10,3000]; % 参数取值上界
   
% The Whale Optimization Algorithm
% function [Leader_score,Leader_pos,Convergence_curve]=WOA(SearchAgents_no,Max_iter,lb,ub,dim,fobj)

% initialize position vector and score for the leader
Leader_pos=zeros(1,dim);
Leader_score=inf; %change this to -inf for maximization problems


%Initialize the positions of search agents
Positions=initialization(SearchAgents_no,dim,ub,lb);

Convergence_curve=zeros(1,Max_iter);

t=0;% Loop counter

% Main loop
while t<Max_iter
    for i=1:size(Positions,1)
        
        % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
       % Position(i,:)=round(Position(i,:));  %%%%%%  自己加上
        % Calculate objective function for each search agent
     %  [u, u_hat, omega] = VMD(signal,  Positions(i,2), tau,  Positions(i,1), DC, init, tol);  %%%%  原始
    [u, u_hat, omega] = VMD(signal, Positions(i,2), tau,  floor(Positions(i,1)), DC, init, tol);
    
     for ii=1:Positions(i,1)
        bao=hilbert(u(ii,:));
        bao=abs(bao);
        p=bao./sum(bao);
        e110(ii,:)=-sum(p.*log10(p));
      end
       fitness=min(e110);
        
       % fitness=fobj(Positions(i,:));
        
        % Update the leader
        if fitness<Leader_score % Change this to > for maximization problem
            Leader_score=fitness; % Update alpha
            Leader_pos=Positions(i,:);
        end
        
    end
    
    a=2-t*((2)/Max_iter); % a decreases linearly fron 2 to 0 in Eq. (2.3)
    
    % a2 linearly dicreases from -1 to -2 to calculate t in Eq. (3.12)
    a2=-1+t*((-1)/Max_iter);
    
    % Update the Position of search agents 
    for i=1:size(Positions,1)
        r1=rand(); % r1 is a random number in [0,1]
        r2=rand(); % r2 is a random number in [0,1]
        
        A=2*a*r1-a;  % Eq. (2.3) in the paper
        C=2*r2;      % Eq. (2.4) in the paper
        
        
        b=1;               %  parameters in Eq. (2.5)
        l=(a2-1)*rand+1;   %  parameters in Eq. (2.5)
        
        p = rand();        % p in Eq. (2.6)
        
        for j=1:size(Positions,2)
            
            if p<0.5   
                if abs(A)>=1
                    rand_leader_index = floor(SearchAgents_no*rand()+1);
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
    t=t+1;
    Convergence_curve(t)=Leader_score;
     [t Leader_score];
end     
bestc=floor(Leader_pos(1,1));  %%%%  自己加上round()
bestg=floor(Leader_pos(1,2));  %%%%  自己加上round()
bestGWOaccuarcy=Leader_score; 

[u, ~, omega] = VMD(signal,  bestg, tau,  bestc, DC, init, tol); 

%out1=DATA;
out1=bestc;   %%%% K
out2=bestg;   %%%%  a
%out4=cz;
out3=u;     %%%% u

% while 1
%     j=0;
% %     [u,~, omega] = VMD(signal, out2, tau, out1, DC, init, tol);
% %     [~, sortIndex] = sort(omega(end,:));
% %     u = u(sortIndex,:);
%     for k=1:bestc
%         a=DFA(u(k,:),s,order,l);
%       if a>0.75
%          j=j+1;
%       end
%     end
%     if j>=J
%        break;
%     else bestc=bestc+1; 
%     end
% end
% x=zeros(1,length(signal));
% for k=1:bestc
%     a=DFA(u(k,:),s,order,l);
%     if a>0.75
%         x=x+u(k,:);
%     end
% end
% DATA=signal-x;
% out4=x;
% out5=DATA;
end




