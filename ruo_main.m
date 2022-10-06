clear;
close all;

t = datetime('now');
save_path = "snr_ser/direct/"+t.Month+"."+t.Day;
if(~exist(save_path,'dir'))
    mkdir(char(save_path));
end
    
channel_choice = 4;
dir_up = "./data_set_final/";
test = [1 0 1 0 1 1 1];
% bias = input('bias(mA): ')

% dp821A = visa('ni','USB0::0x1AB1::0x0E11::DP8E163250125::INSTR');
% % dp821A = visa('ni','ASRL3:INSTR');
% fopen(dp821A);
% % fprintf(dp821A,'*IDN?');
% % aa = fscanf(dp821A);
% % display(aa);
% fprintf(dp821A,':OUTP:OCP:VAL CH2,1.500');
% fprintf(dp821A,':OUTP:OCP:VAL? CH2');
% aa = fscanf(dp821A);
% display(aa);
% fprintf(dp821A,':OUTP:OCP CH2,ON');
% fprintf(dp821A,':OUTP CH2,ON');
% fprintf(dp821A,':APPL CH1,24,1.000');
% fprintf(dp821A,':OUTP CH1,ON');

bias_name = 780;
pause(0.5);

data_length = 10000;
ls_order = 50;

origin_rate = 25e6; 
bw = origin_rate/2;         % baseband bandwidth
f_rate = 160e6;
d_rate = 150e6;

times = d_rate/origin_rate;
num_of_windows = 100;

% pilot = pilot_gen([0 0 0 1 1 1 1]);
% pilot = pilot_gen([1 1 1 0 0 0 1 1 1]);
pilot = ruo_pilot_gen([1 1 1 1 0 0 1 1 0 1 1]);
pilot_ps = bandpower(pilot);
pilot_length = length(pilot);
pilot_bpsk = pilot*2-1;

filter_ord = 1000;     % filter order used in function sam_rate_con
rp = 0.00057565;      
rst = 1e-4;       % filter parameter used in function sam_rate_con

llcm_transmit = lcm(origin_rate,f_rate);
upf_transmit = llcm_transmit/origin_rate;  % upsampling parameters
dof_transmit = llcm_transmit/f_rate;
ups_rate_transmit = origin_rate*upf_transmit;
filter_transmit = firceqrip(filter_ord,bw/(ups_rate_transmit/2),[rp rst],'passedge'); % Impulse function
filter_transmit = filter_transmit./norm(filter_transmit,2)*sqrt(bw/ups_rate_transmit*2); % filter used in function sam_rate_con in transmit part

llcm_receive = lcm(f_rate,d_rate);
upf_receive = llcm_receive/f_rate;  % upsampling parameters
dof_receive = llcm_receive/d_rate;
ups_rate_receive = f_rate*upf_receive;
filter_receive = firceqrip(filter_ord,bw/(ups_rate_receive/2),[rp rst],'passedge'); % Impulse function
filter_receive = filter_receive./norm(filter_receive,2)*sqrt(bw/ups_rate_receive*2); % filter used in function sam_rate_con in receive part

% origin_rate = gpuArray(double(origin_rate));
% f_rate = gpuArray(double(f_rate));
% d_rate = gpuArray(double(d_rate));
% num_of_windows = gpuArray(double(num_of_windows));
% times = gpuArray(double(times));
% pilot_length = gpuArray(double(pilot_length));
% delay = gpuArray(double(delay));
% filter_transmit = gpuArray(double(filter_transmit));
% filter_receive = gpuArray(double(filter_receive));
% upf_transmit = gpuArray(double(upf_transmit));
% dof_transmit = gpuArray(double(dof_transmit));
% upf_receive = gpuArray(double(upf_receive));
% dof_receive = gpuArray(double(dof_receive));
% h_channel = gpuArray(double(h_channel));
% h_channel_delay = gpuArray(double(h_channel_delay));

snr_begin = 3;
snr_end = 22;
fprintf('add zero,ls order=%d,pilot length=%d .\n',ls_order,pilot_length);
for snr = 22:snr_end
    looptime = 0;
    snr_sum = 0;
    errornum_zf = 0;
    errornum_mmse = 0;
    errornum_ls = 0;
    fprintf('snr = %d .\n', snr);

    while(errornum_ls <= 30 || looptime < 2000)
%     while(errornum_zf <= 100 || errornum_mmse <= 100 || looptime < 50)
        if mod(looptime,1000) == 0
            ruo_zero_send;
            pn = bandpower(noise)*(8192^2);
        end
        ps = pn*10^(snr/10);
        looptime = looptime+1;
        
        ruo_pam4_send;
          
        signal_received = ruo_sam_rate_con(signal_pass_channel,filter_receive,upf_receive,dof_receive);

        [fin_syn_point,coar_syn_point] = ruo_signal_syn(origin_rate,d_rate,signal_ori,signal_received,num_of_windows);
            
        signal_demod_ls = ruo_signal_equal(signal_ori,signal_received,times,coar_syn_point,pilot_length,num_of_windows,{'ls',ls_order,pilot_ps,data_mpam_ps});

%         signal_demod_ls = gather(signal_demod_ls);

        data_demod_ls = signal_demod_ls(pilot_length+11:end);

        if length(data_demod_ls) < length(data)
            errornum_ls_loop = sum(data_demod_ls ~= data(1:length(data_demod_ls)));
        else
            errornum_ls_loop = sum(data_demod_ls(1:length(data)) ~= data);
        end
        
        errornum_ls = errornum_ls+errornum_ls_loop;

        ser_ls = errornum_ls/(looptime*length(data_demod_ls));

        if mod(looptime,2) == 0
           fprintf(' %f times, snr = %d , data num = %d ,ls error num = %d .\n',looptime,snr,length(data_demod_ls),errornum_ls_loop);
           fprintf(' %f times, snr = %d , total ls error num = %d,ls error rate = %.6g .\n',looptime,snr,errornum_ls,ser_ls);

%            fprintf(' %f times, snr = %f , zf error num = %f , mmse error num = %f .\n',looptime,snr,errornum_zf,errornum_mmse);
%            fprintf(' %f times, snr = %f , zf error rate = %f , mmse error rate = %f .\n',looptime,snr,ser_zf,ser_mmse);
        end
    end
    
%     ser_ls = gather(ser_ls);

    if snr == snr_begin
        fsnr = fopen(save_path+"/snr.txt",'w');
        fser = fopen(save_path+"/ser.txt",'w');
        fprintf(fsnr,'add zero ,pilot length  = %.8f , ls order  = %.8f \r\n',pilot_length,ls_order);
        fprintf(fser,'add zero ,pilot length  = %.8f , ls order  = %.8f \r\n',pilot_length,ls_order);
    else
        fsnr = fopen(save_path+"snr.txt",'a');
        fser = fopen(save_path+"ser.txt",'a');
    end
    fprintf(fsnr,'%.8f \r\n',snr);
    fprintf(fser,'%.6g \r\n',ser_ls);
    fclose(fsnr);
    fclose(fser);

end


