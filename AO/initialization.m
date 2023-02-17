%% �Ӻ������ڳ�ʼ����Ⱥλ��
function X=initialization(N,Dim,UB,LB)

B_no= size(UB,2); % ����ά��

% ��ά��
if B_no==1
    % ��ʼ��λ��
    X=rand(N,Dim).*(UB-LB)+LB;
end

% ��ά��
if B_no>1
    for i=1:Dim
        Ub_i=UB(i);
        Lb_i=LB(i);
        % ��ʼ��λ��
        X(:,i)=rand(N,1).*(Ub_i-Lb_i)+Lb_i;
    end
end