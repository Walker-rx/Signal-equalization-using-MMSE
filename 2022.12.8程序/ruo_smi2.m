clear;
close all;

t = datetime('now');
save_path = "snr_ser/direct/"+t.Month+"."+t.Day+"/"+t.Hour+"_smi";
if(~exist(save_path,'dir'))
    mkdir(char(save_path));
end
    
channel_choice = 4;
dir_up = "./data_set_final/";
test = [1 0 1 0 1 1 1];

bias_name = 780;
pause(0.5);

data_length = 100000;
zero_length = 10000;
ls_order = 50;

origin_rate = 25e6; 
bw = origin_rate/2;         % baseband bandwidth
f_rate = 160e6;
d_rate = 150e6;

times = d_rate/origin_rate;
num_of_windows = 100;

% pilot = pilot_gen([0 0 0 1 1 1 1]);
% pilot = pilot_gen([1 1 1 0 0 0 1 1 1]);
pilot = ruo_pilot_gen([1 1 1 1 1 1 0 1 0 1 0 1 1]);
% pilot = ruo_pilot_gen([1 1 1 1 1 1 0 1 0 0 1]);
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

filter_channel = firceqrip(filter_ord,bw*5/(ups_rate_transmit/2),[rp rst],'passedge'); % Impulse function
filter_channel = filter_channel./norm(filter_channel,2)*sqrt(bw*5/ups_rate_transmit*2); % filter used in function sam_rate_con in channel part

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

amp_begin = 20;
amp_end = 50;
fprintf('add zero,ls order=%d,pilot length=%d .\n',ls_order,pilot_length);
for amp = 35:amp_end
    looptime = 0;
    snr = amp-15;
    ps_smi = 0;
    pn_smi = 0;
    errornum_ls = 0;
    total_length = 0;
    fprintf('amp = %d .\n', amp);
    
    while(errornum_ls <= 30 || looptime < 2000)
  
        looptime = looptime+1;
        
        ruo_pam4_send;
        
        signal_received = ruo_sam_rate_con(signal_pass_channel,filter_receive,upf_receive,dof_receive);

        [fin_syn_point,coar_syn_point] = ruo_signal_syn(origin_rate,d_rate,signal_ori,signal_received,num_of_windows);
        
        signal_downsample = signal_received(fin_syn_point:times:end);
        noise = signal_downsample(pilot_length+100:pilot_length+6000);
        pn_channel = bandpower(noise)*8192^2;
        ps_channel = pn_channel*10^(snr/10);
        
        channel_smi_bw = ruo_channel_coefficient(signal_ori,signal_received,times,coar_syn_point,pilot_length,num_of_windows);
        channel_smi_up = ruo_sam_rate_con(channel_smi_bw,filter_channel,upf_transmit,dof_transmit);
        [~,index_bw] = max(channel_smi_bw);
        [~,index_up] = max(channel_smi_up);
        index_correct = ceil(index_bw*f_rate/origin_rate);
        index_bias = index_correct - index_up;
        channel_smi = [ zeros(1,index_bias) , channel_smi_up ];
        channel_smi = channel_smi./norm(channel_smi,2)*sqrt(length(channel_smi));
        channel_delay_smi = length(channel_smi);
        
        signal_upsample_smi = signal_upsample;
        signal_upsample_smi = signal_upsample_smi./norm(signal_upsample_smi,2)*sqrt(length(signal_upsample_smi))*sqrt(ps_channel);
        signal_send_smi = [rand([1,512]) zeros(1,10) signal_upsample_smi];
        
        signal_passtap_smi = conv(signal_send_smi,channel_smi);
        signal_passtap_smi = signal_passtap_smi((channel_delay_smi+1)/2:length(signal_passtap_smi)-(channel_delay_smi-1)/2);
        signal_passtap_smi = signal_passtap_smi./norm(signal_passtap_smi,2)*sqrt(length(signal_passtap_smi))*sqrt(ps_channel);
        signal_passchannel_smi = awgn(signal_passtap_smi,snr,'measured');
        signal_passchannel_smi = signal_passchannel_smi/8192;
        
        signal_received_smi = ruo_sam_rate_con(signal_passchannel_smi,filter_receive,upf_receive,dof_receive);
       
        [fin_syn_point_smi,coar_syn_point_smi] = ruo_signal_syn(origin_rate,d_rate,signal_ori,signal_received_smi,num_of_windows);
        
        signal_downsample_smi = signal_received_smi(fin_syn_point_smi:times:end);
        noise_smi = signal_downsample_smi(pilot_length+100:pilot_length+6000);
        pn_loop_smi = bandpower(noise_smi);
        data_received_smi = signal_downsample_smi(pilot_length+zero_length+1:pilot_length+zero_length+1+data_length-1);
        p_smi = bandpower(data_received_smi);
        ps_loop_smi = p_smi - pn_loop_smi;
        
        signal_demod_ls_smi = ruo_signal_equal(signal_ori,signal_received_smi,times,coar_syn_point_smi,pilot_length,zero_length,num_of_windows,{'ls',ls_order});

%         signal_demod_ls = gather(signal_demod_ls);

        data_demod_ls_smi = signal_demod_ls_smi(pilot_length+zero_length+1:end);

        if length(data_demod_ls_smi) < length(data)
            compare_length = length(data_demod_ls_smi);
            total_length = total_length + compare_length;
            errornum_ls_loop = sum(data_demod_ls_smi ~= data(1:compare_length));
            error_location = find(data_demod_ls_smi ~= data(1:compare_length));
        else
            compare_length = length(data);
            total_length = total_length + compare_length;
            errornum_ls_loop = sum(data_demod_ls_smi(1:compare_length) ~= data);
            error_location = find(data_demod_ls_smi(1:compare_length) ~= data);
        end
        
        errornum_ls = errornum_ls+errornum_ls_loop;      
        ser_ls = errornum_ls/total_length;

        ps_smi = ps_smi + ps_loop_smi;
        pn_smi = pn_smi + pn_loop_smi;
        snr_ls_smi = 10*log10(ps_smi/pn_smi);
        if mod(looptime,2) == 0
           fprintf('smi2 ,  %f times, amp = %d , data num = %d ,ls error num = %d .\n',looptime,amp,length(data_demod_ls_smi),errornum_ls_loop);
           fprintf('smi2 ,  %f times, snr = %d , total ls error num = %d,ls error rate = %.6g .\n',looptime,snr_ls_smi,errornum_ls,ser_ls);
%            disp(["error location = ",error_location]);
%            disp(["correct = ",data(error_location)]);
%            disp(["false = ",data_demod_ls_smi(error_location)]);
        end
          
    end
    
%     ser_ls = gather(ser_ls);

    if amp == amp_begin
        fsnr = fopen(save_path+"/snr_smi2.txt",'w');
        fser = fopen(save_path+"/ser_smi2.txt",'w');
        fprintf(fsnr,'add zero ,simulation ,pilot length  = %.8f , ls order  = %.8f \r\n',pilot_length,ls_order);
        fprintf(fser,'add zero ,simulation ,pilot length  = %.8f , ls order  = %.8f \r\n',pilot_length,ls_order);
    else
        fsnr = fopen(save_path+"/snr_smi2.txt",'a');
        fser = fopen(save_path+"/ser_smi2.txt",'a');
    end
    fprintf(fsnr,'%.8f \r\n',snr_ls_smi);
    fprintf(fser,'%.6g \r\n',ser_ls);
    fclose(fsnr);
    fclose(fser);

end


