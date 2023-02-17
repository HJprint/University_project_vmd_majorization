function ret=Decode(lenchrom,bound,code,opts)
% ��������Ⱦɫ����н���
% lenchrom   input : Ⱦɫ�峤��
% bound      input : ����ȡֵ��Χ
% code       input ������ֵ
% opts       input : ���뷽����ǩ
% ret        output: Ⱦɫ��Ľ���ֵ
switch opts
    case 'binary' % binary coding
        for i=length(lenchrom):-1:1
        data(i)=bitand(code,2^lenchrom(i)-1);  %����ʮλ��Ȼ�󽫵�ʮλת����ʮ����������data(i)����
        code=(code-data(i))/(2^lenchrom(i));   %��ʮλ���㣬Ȼ������ʮλ
        end
        ret=bound(:,1)'+data./(2.^lenchrom-1).*(bound(:,2)-bound(:,1))';  %�ֶν��룬��ʵ����������ʽ����ret��
        
    case 'grey'   % grey coding
        for i=sum(lenchrom):-1:2
            code=bitset(code,i-1,bitxor(bitget(code,i),bitget(code,i-1)));
        end
        for i=length(lenchrom):-1:1
        data(i)=bitand(code,2^lenchrom(i)-1);
        code=(code-data(i))/(2^lenchrom(i));
        end
        ret=bound(:,1)'+data./(2.^lenchrom-1).*(bound(:,2)-bound(:,1))'; %�ֶν��룬��ʵ����������ʽ����ret��
        
    case 'float'  % float coding
        ret=code; %���������Ǳ�������ʵ��������������ret��
end