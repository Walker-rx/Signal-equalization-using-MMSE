clear
% origin_rate = 25e6; 
origin_rate = 1.8e6;
bw = origin_rate/2;         % baseband bandwidth
f_rate = 160e6;
% d_rate = 150e6;
d_rate = origin_rate*6;

err = 1e20;
for n = 1:100
    llcm_transmit_tmp=(f_rate*n)/round(f_rate*n/origin_rate);
    err_tmp = abs(llcm_transmit_tmp - origin_rate);
    if (err_tmp < err) && ((f_rate*n)/llcm_transmit_tmp <= 200)
        err = err_tmp;
        N = n
        new_orignal_rate = llcm_transmit_tmp;
        (f_rate*n)/llcm_transmit_tmp
    end

end
