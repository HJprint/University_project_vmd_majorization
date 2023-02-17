%% 子函数用于初始化种群位置
function X=initialization(N,Dim,UB,LB)

B_no= size(UB,2); % 搜索维度

% 单维度
if B_no==1
    % 初始化位置
    X=rand(N,Dim).*(UB-LB)+LB;
end

% 多维度
if B_no>1
    for i=1:Dim
        Ub_i=UB(i);
        Lb_i=LB(i);
        % 初始化位置
        X(:,i)=rand(N,1).*(Ub_i-Lb_i)+Lb_i;
    end
end