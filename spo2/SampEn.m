
function SampEnVal = SampEn(data, m, r)
%SAMPEN  计算时间序列data的样本熵
%        data为输入数据序列
%        m为初始分段，每段的数据长度
%        r为阈值
% $Author: lskyp
% $Date:   2010.6.20
% Orig Version: V1.0--------分开计算长度为m的序列和长度为m+1的序列
%                           这一版的计算有些问题，需要注意两个序列总数都要为N-m
% Modi Version: V1.1--------综合计算，计算距离时通过矩阵减法完成，避免重循环
% V1.1 Modified date: 2010.6.23
data = data(:)';
N = length(data);
Nkx1 = 0;
Nkx2 = 0;
% 分段计算距离，x1为长度为m的序列，x2为长度为m+1的序列
for k = N - m:-1:1
    x1(k, :) = data(k:k + m - 1);
    x2(k, :) = data(k:k + m);
end
for k = N - m:-1:1
    % x1序列计算
    % 统计距离，由于每行都要与其他行做减法，因此可以先将该行复制为N-m的矩阵，然后
    % 与原始x1矩阵做减法，可以避免两重循环，增加效率
    x1temprow = x1(k, :);
    x1temp    = ones(N - m, 1)*x1temprow;
    % 可以使用repmat函数完成上面的语句，即
    % x1temp = repmat(x1temprow, N - m, 1);
    % 但是效率不如上面的矩阵乘法
    % 计算距离，每一行元素相减的最大值为距离
    dx1(k, :) = max(abs(x1temp - x1), [], 2)';
    % 模板匹配数
    Nkx1 = Nkx1 + (sum(dx1(k, :) < r) - 1)/(N - m - 1);
    
    % x2序列计算，和x1同样方法
    x2temprow = x2(k, :);
    x2temp    = ones(N - m, 1)*x2temprow;
    dx2(k, :) = max(abs(x2temp - x2), [], 2)';
    Nkx2      = Nkx2 + (sum(dx2(k, :) < r) - 1)/(N - m - 1);
end
% 平均值
Bmx1 = Nkx1/(N - m);
Bmx2 = Nkx2/(N - m);
% 样本熵
SampEnVal = -log(Bmx2/Bmx1);
