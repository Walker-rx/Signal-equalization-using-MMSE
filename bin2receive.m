function Rx=bin2receive(ch_num)
Fid = fopen("./binUDP/receive.bin",'r');
Resig = fread(Fid,'int16');
fclose('all');
Resig = reshape(Resig,4,length(Resig)/4);
Rx = Resig(ch_num,:);
Rx = Rx';
