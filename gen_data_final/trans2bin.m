function trans2bin(X)

if(~exist("./binUDP/",'dir'))
        mkdir(char("./binUDP/"));
end
Fid = fopen("./binUDP/transmit.bin",'w');
fwrite(Fid,X,'int16');
fclose('all');
stop_addr = length(X)/704*176;
Fid2 = fopen("./binUDP/ddr_stop_addr.bin",'w');
fwrite(Fid2,stop_addr,'int32');
fclose('all');