%% 子函数 位置越界限制
function s=simplebounds(s,Lb,Ub)
%下限
ns_tmp=s;
Lb = Lb.*ones(1, length(s));
I=ns_tmp<Lb;
ns_tmp(I)=Lb(I);

%上限
Ub = Ub.*ones(1, length(s));
J=ns_tmp>Ub;
ns_tmp(J)=Ub(J);

s=ns_tmp;
end