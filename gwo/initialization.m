%___________________________________________________________________%
%  Grey Wolf Optimizer (GWO) source codes version 1.0               %
%                                                                   %
%  Developed in MATLAB R2011b(7.13)                                 %
%                                                                   %
%  Author and programmer: Seyedali Mirjalili                        %
%                                                                   %
%         e-Mail: ali.mirjalili@gmail.com                           %
%                 seyedali.mirjalili@griffithuni.edu.au             %
%                                                                   %
%       Homepage: http://www.alimirjalili.com                       %
%                                                                   %
%   Main paper: S. Mirjalili, S. M. Mirjalili, A. Lewis             %
%               Grey Wolf Optimizer, Advances in Engineering        %
%               Software , in press,                                %
%               DOI: 10.1016/j.advengsoft.2013.12.007               %
%                                                                   %
%___________________________________________________________________%

% This function initialize the first population of search agents
function Positions=initialization(SearchAgents_no,dim,ub,lb)%SearchAgents_no=10
%SearchAgents_no=10; dim=2; ub=[10,5000]; lb=[2,100];
Boundary_no= size(ub,2); %ub的第二个维度的长度。 numnber of boundaries

%如果所有变量的边界相等，用户输入一个 If the boundaries of all variables are equal and user enter a signle
% number for both ub and lb
if Boundary_no==1
    Positions=rand(SearchAgents_no,dim).*(ub-lb)+lb;%一个10*2的随机限制矩阵
end

% 如果所有变量的边界不相等If each variable has a different lb and ub
%因为参数取值上界和下界有两个范围2-10和100-5000，所以对于不同的限制有不同的位置
if Boundary_no>1
    for i=1:dim%dim=2
        ub_i=ub(i);
        lb_i=lb(i);
        Positions(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;%计算两个界限的初始化位置
        %Positions为一个10*2的随机限制矩阵
    end
end