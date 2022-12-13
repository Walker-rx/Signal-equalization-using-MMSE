function X = ruo_TxDataSort(Tx)
% Tx: cell数组
ch_num = length(Tx);
odd_zeros = 1407; %奇数补零
even_zeros = 1408; %偶数补零
if ch_num==1
    Xa = Tx{1};
    Xb = Tx{1};
    Xc = Tx{1};
    Xd = Tx{1};
elseif ch_num ==2
    Xa = Tx{1};% DA2 in hardware
    Xb = Tx{2}; % DA3 in hardware
    Xc = Tx{1}; % DA4 in hardware
    Xd = Tx{2}; % DA1 in hardwareS
elseif ch_num ==3
    Xa = Tx{1};
    Xb = Tx{2};
    Xc = Tx{3};
    Xd = Tx{1};
elseif ch_num ==4
    Xa = Tx{1};
    Xb = Tx{2};
    Xc = Tx{3};
    Xd = Tx{4};
elseif ch_num ==0
    error("输入参数不足")
else
    error("输入参数太多")
end
Xa = reshape(Xa,[1,length(Xa)]);
Xb = reshape(Xb,[1,length(Xb)]);
Xc = reshape(Xc,[1,length(Xc)]);
Xd = reshape(Xd,[1,length(Xd)]);
if mod(length(Xa),2)==1
    Xa = [Xa,zeros(1,odd_zeros)];
    Xb = [Xb,zeros(1,odd_zeros)];
    Xc = [Xc,zeros(1,odd_zeros)];
    Xd = [Xd,zeros(1,odd_zeros)];
else
    Xa = [Xa,zeros(1,even_zeros)];
    Xb = [Xb,zeros(1,even_zeros)];
    Xc = [Xc,zeros(1,even_zeros)];
    Xd = [Xd,zeros(1,even_zeros)];
end 
N = 4*length(Xa);
tmp_matr = zeros(N/8,8);
X = zeros(N,1);
for i = 1:N/8
    for k = 1:8
        tmp_matr(i,k) = (k==1)*Xa(2*i) + (k==2)*Xa(2*i-1) + (k==3)*Xb(2*i) + (k==4)*Xb(2*i-1) + (k==5)*Xc(2*i) + (k==6)*Xc(2*i-1)...
            + (k==7)*Xd(2*i) + (k==8)*Xd(2*i-1);
    end
end
for i = 1:N/8
    for k = 1:8
        X((i-1)*8+k) = tmp_matr(i,k);
    end
end
%尾部补零防止丢失尾部数据
% X = [X;zeros(2000,1)];