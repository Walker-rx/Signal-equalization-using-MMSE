function Rx=bin2receive_test(ch_num)
%%
 % ch_num 1 (AD6 in hardware)
 % ch_num 2  (AD8 in hardware) 
 % ch_num 3 (AD2 in hardware)
 % ch_num 4 (AD4 in hardware)
%%
Fid = fopen("./data_bin/receive_test.bin",'r');
Resig = fread(Fid,'int16');
fclose('all');
Resig = reshape(Resig,4,length(Resig)/4);
Rx = Resig(ch_num,:);
Rx = Rx';
