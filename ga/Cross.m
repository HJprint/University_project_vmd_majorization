function ret=Cross(pcross,lenchrom,chrom,sizepop,bound)
%��������ɽ������
% pcorss                input  : �������
% lenchrom              input  : Ⱦɫ��ĳ���
% chrom     input  : Ⱦɫ��Ⱥ
% sizepop               input  : ��Ⱥ��ģ
% ret                   output : ������Ⱦɫ��
 for i=1:sizepop  %ÿһ��forѭ���У����ܻ����һ�ν��������Ⱦɫ�������ѡ��ģ�����λ��Ҳ�����ѡ��ģ�%������forѭ�����Ƿ���н���������ɽ�����ʾ�����continue���ƣ�
     % ���ѡ������Ⱦɫ����н���
     pick=rand(1,2);
     while prod(pick)==0
         pick=rand(1,2);
     end
     index=ceil(pick.*sizepop);
     % ������ʾ����Ƿ���н���
     pick=rand;
     while pick==0
         pick=rand;
     end
     if pick>pcross
         continue;
     end
     flag=0;
     while flag==0
         % ���ѡ�񽻲�λ
         pick=rand;
         while pick==0
             pick=rand;
         end
         pos=ceil(pick.*sum(lenchrom)); %���ѡ����н����λ�ã���ѡ��ڼ����������н��棬ע�⣺����Ⱦɫ�彻���λ����ͬ
         pick=rand; %���濪ʼ
         v1=chrom(index(1),pos);
         v2=chrom(index(2),pos);
         chrom(index(1),pos)=pick*v2+(1-pick)*v1;%  pick*v2+(1-pick)*v1=v2 �����Ƴ�v1=v2  �൱����
         chrom(index(2),pos)=pick*v1+(1-pick)*v2; %�������
         flag1=test(lenchrom,bound,chrom(index(1),:));  %����Ⱦɫ��1�Ŀ�����
         flag2=test(lenchrom,bound,chrom(index(2),:));  %����Ⱦɫ��2�Ŀ�����
         if   flag1*flag2==0
             flag=0;
         else flag=1;
         end    %�������Ⱦɫ�岻�Ƕ����У������½���
     end
 end
ret=chrom;