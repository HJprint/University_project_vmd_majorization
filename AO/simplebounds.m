%% �Ӻ��� λ��Խ������
function s=simplebounds(s,Lb,Ub)
%����
ns_tmp=s;
Lb = Lb.*ones(1, length(s));
I=ns_tmp<Lb;
ns_tmp(I)=Lb(I);

%����
Ub = Ub.*ones(1, length(s));
J=ns_tmp>Ub;
ns_tmp(J)=Ub(J);

s=ns_tmp;
end