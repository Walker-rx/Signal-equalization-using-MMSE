clear;
close all;

t = datetime('now');
save_path = "snr_ser/direct/"+t.Month+"."+t.Day+'/delete2';
if(~exist(save_path,'dir'))
    mkdir(char(save_path));
end

channel_choice = 1; % % 1 is direct channel
channel_choice_inf = 4;
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

M = 4;
data_length = 10000;
zero_length = 3000;
ls_order = 50;

origin_rate = 25e6; 
bw = origin_rate/2;         % baseband bandwidth
f_rate = 160e6;
d_rate = 150e6;

times = d_rate/origin_rate;
num_of_windows = 100;

% pilot = pilot_gen([0 0 0 1 1 1 1]);
% pilot = pilot_gen([1 1 1 0 0 0 1 1 1]);
pilot = ruo_pilot_gen([1 1 1 1 1 1 0 1 0 0 1]);    %  2^11 = 2048
% pilot = ruo_pilot_gen([1 1 1 1 1 1 0 1 0 1 0 1 1]);   %  2^13 = 8192
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
rng(1);
amp_begin = 50;
amp_end = 53;
fprintf('add zero,ls order=%d,pilot length=%d .\n',ls_order,pilot_length);
for amp = amp_begin:amp_end
    looptime = 0;
    ps = 0;
    pn = 0;
    ps_beforecorrect = 0;
    pn_beforecorrect = 0;
    errornum_zf = 0;
    errornum_mmse = 0;
    errornum_ls = 0;
    errornum_ls_beforecorrect = 0;
    total_length_nocorrect = 0;
    total_length_correct = 0;
    total_length = 0;
    txerror_num = 0;
    fprintf('amp = %d .\n', amp);
    
    while(errornum_ls <= 30 || looptime < 1000)
%     while(errornum_zf <= 100 || errornum_mmse <= 100 || looptime < 50)
        
        looptime = looptime+1;
        
        data = randi([0,M-1],[1,data_length]);
        data_mpam = real(pammod(data,M));
        data_mpam_ps = bandpower(data_mpam);
        
        unmod = [pilot,zeros(1,zero_length),data];
        signal_ori = [pilot_bpsk,zeros(1,zero_length),data_mpam];
        signal_upsample = ruo_sam_rate_con(signal_ori,filter_transmit,upf_transmit,dof_transmit);       
        signal_send_tmp = signal_upsample./norm(signal_upsample,2)*sqrt(length(signal_upsample))*100*1.1^(amp);
        signal_send_tmp_inf = signal_upsample./norm(signal_upsample,2)*sqrt(length(signal_upsample))*100*1.1^(amp_end);
        signal_rand = rand([1,512]);
        signal_send = [signal_rand zeros(1,10) signal_send_tmp];
        signal_send_inf = [signal_rand zeros(1,10) signal_send_tmp_inf];
      
        ruo_pam4_send;
        ruo_pam4_send_correct;
        
        %% Locating the transmission error point
        signal_received_inf = ruo_sam_rate_con(signal_pass_channel_inf,filter_receive,upf_receive,dof_receive);
        
        [fin_syn_point_inf,coar_syn_point_inf] = ruo_signal_syn(origin_rate,d_rate,signal_ori,signal_received_inf,num_of_windows);
        
        signal_demod_ls_inf = ruo_signal_equal(signal_ori,signal_received_inf,times,coar_syn_point_inf,pilot_length,zero_length,num_of_windows,{'ls',ls_order});
        
%         data_demod_ls_inf = signal_demod_ls_inf(pilot_length+zero_length+1:end);
        
        if length(signal_demod_ls_inf) < length(unmod)
            compare_length_inf = length(signal_demod_ls_inf);
            error_location_inf = find(signal_demod_ls_inf ~= unmod(1:compare_length_inf));
        else
            compare_length_inf = length(unmod);
            error_location_inf = find(signal_demod_ls_inf(1:compare_length_inf) ~= unmod);
        end
        txerror_num_loop = length(error_location_inf) ;
        txerror_num = txerror_num + txerror_num_loop ;
        txerror_num_average = txerror_num/looptime ;
        %% Calculating the SER without correction
        signal_received = ruo_sam_rate_con(signal_pass_channel,filter_receive,upf_receive,dof_receive);       
        [fin_syn_point,coar_syn_point] = ruo_signal_syn(origin_rate,d_rate,signal_ori,signal_received,num_of_windows);
        
        signal_received_correct = ruo_sam_rate_con(signal_pass_channel_correct,filter_receive,upf_receive,dof_receive);
        [fin_syn_point_correct,coar_syn_point_correct] = ruo_signal_syn(origin_rate,d_rate,signal_ori,signal_received_correct,num_of_windows);
        
%         signal_downsample_nocorrect = signal_received(fin_syn_point:times:end);
%         noise_nocorrect = signal_downsample_nocorrect(pilot_length+100:pilot_length+100+zero_length/2-1);
%         pn_loop_nocorrect = bandpower(noise_nocorrect);
%         data_received_nocorrect = signal_downsample_nocorrect(pilot_length+zero_length+1:pilot_length+zero_length+1+data_length-1);
%         p_nocorrect = bandpower(data_received_nocorrect);
%         ps_loop_nocorrect = p_nocorrect - pn_loop_nocorrect;
%         
%         signal_demod_ls_nocorrect = ruo_signal_equal(signal_ori,signal_received,times,coar_syn_point,pilot_length,zero_length,num_of_windows,{'ls',ls_order});
%         data_demod_ls_nocorrect = signal_demod_ls_nocorrect(pilot_length+zero_length+1:end);
%         
%         if length(data_demod_ls_nocorrect) < data_length
%             compare_length = length(data_demod_ls_nocorrect);
%             total_length_nocorrect = total_length_nocorrect + compare_length;
%             errornum_ls_loop_nocorrect = sum(data_demod_ls_nocorrect ~= data(1:compare_length));
%             error_location_nocorrect = find(data_demod_ls_nocorrect ~= data(1:compare_length));
%         else
%             compare_length = data_length;
%             total_length_nocorrect = total_length_nocorrect + compare_length;
%             errornum_ls_loop_nocorrect = sum(data_demod_ls_nocorrect(1:compare_length) ~= data);
%             error_location_nocorrect = find(data_demod_ls_nocorrect(1:compare_length) ~= data);
%         end
%         
%         signal_downsample_correct = signal_received_correct(fin_syn_point_correct:times:end);
%         noise_correct = signal_downsample_correct(pilot_length+100:pilot_length+100+zero_length/2-1);
%         pn_loop_correct = bandpower(noise_correct);
%         data_received_correct = signal_downsample_correct(pilot_length+zero_length+1:pilot_length+zero_length+1+data_length-1);
%         p_correct = bandpower(data_received_correct);
%         ps_loop_correct = p_correct - pn_loop_correct;
%         
%         signal_demod_ls_correct = ruo_signal_equal(signal_ori,signal_received_correct,times,coar_syn_point_correct,pilot_length,zero_length,num_of_windows,{'ls',ls_order});
%         data_demod_ls_correct = signal_demod_ls_correct(pilot_length+zero_length+1:end);
%         
%         if length(data_demod_ls_correct) < data_length
%             compare_length = length(data_demod_ls_correct);
%             total_length_correct = total_length_correct + compare_length;
%             errornum_ls_loop_correct = sum(data_demod_ls_correct ~= data(1:compare_length));
%             error_location_correct = find(data_demod_ls_correct ~= data(1:compare_length));
%         else
%             compare_length = data_length;
%             total_length_correct = total_length_correct + compare_length;
%             errornum_ls_loop_correct = sum(data_demod_ls_correct(1:compare_length) ~= data);
%             error_location_correct = find(data_demod_ls_correct(1:compare_length) ~= data);
%         end
%         
%         total_length_beforecorrect = total_length_correct + total_length_nocorrect;
%         errornum_ls_beforecorrect = errornum_ls_beforecorrect + errornum_ls_loop_correct + errornum_ls_loop_nocorrect;
%         ps_beforecorrect = ps_beforecorrect + ps_loop_correct + ps_loop_nocorrect;
%         pn_beforecorrect = pn_beforecorrect + pn_loop_correct + pn_loop_nocorrect;
%         
%         ser_ls_beforecorrect = errornum_ls_beforecorrect/total_length_beforecorrect;
%         snr_ls_beforecorrect = 10*log10(ps_beforecorrect/pn_beforecorrect);
        %% Replacing the transmission error point
        signal_received_tmp = signal_received;
        replace_loc = error_location_inf*times + fin_syn_point - 1;
        replace_loc_correct = error_location_inf*times +fin_syn_point_correct - 1;
        replace_length = 8*times;
        for i = 1:txerror_num_loop
            signal_received(replace_loc(i)-replace_length : replace_loc(i)+replace_length) = signal_received_correct(replace_loc_correct(i)-replace_length : replace_loc_correct(i)+replace_length);
        end       
        %% Calculating the SER after correction
        signal_downsample = signal_received(fin_syn_point:times:end);
        noise = signal_downsample(pilot_length+100:pilot_length+100+zero_length/2-1);
        pn_loop = bandpower(noise);
        data_received = signal_downsample(pilot_length+zero_length+1:pilot_length+zero_length+1+data_length-1);
        p = bandpower(data_received);
        ps_loop = p - pn_loop;
             
        signal_demod_ls = ruo_signal_equal(signal_ori,signal_received,times,coar_syn_point,pilot_length,zero_length,num_of_windows,{'ls',ls_order});

%         signal_demod_ls = gather(signal_demod_ls);

        data_demod_ls = signal_demod_ls(pilot_length+zero_length+1:end);

        if length(data_demod_ls) < data_length
            compare_length = length(data_demod_ls);
            total_length = total_length + compare_length;
            errornum_ls_loop = sum(data_demod_ls ~= data(1:compare_length));
            error_location = find(data_demod_ls ~= data(1:compare_length));
        else
            compare_length = data_length;
            total_length = total_length + compare_length;
            errornum_ls_loop = sum(data_demod_ls(1:compare_length) ~= data);
            error_location = find(data_demod_ls(1:compare_length) ~= data);
        end
        
        errornum_ls = errornum_ls+errornum_ls_loop;
        ps = ps + ps_loop;
        pn = pn + pn_loop;
        
        ser_ls = errornum_ls/total_length;
        snr_ls = 10*log10(ps/pn);
        
        snr_ls_beforecorrect = 0;
        ser_ls_beforecorrect = 0; 
        if mod(looptime,4) == 0
           fprintf(' 4pam , amp = %d , %d times, data num = %d ,average tx error num = %d .\n',amp,looptime,length(data_demod_ls),txerror_num_average);
           fprintf(' %f times, before, snr = %d , total ls error num = %d , ls error rate = %.6g .\n',looptime,snr_ls_beforecorrect,errornum_ls_beforecorrect/2,ser_ls_beforecorrect);
           fprintf(' %f times, after, snr = %d , loop ls error num = %d , total ls error num = %d , ls error rate = %.6g .\n',looptime,snr_ls,errornum_ls_loop,errornum_ls,ser_ls);
%            disp(["error location = ",error_location]);
%            disp(["correct = ",data(error_location)]);
%            disp(["false = ",data_demod_ls(error_location)]);
%            fprintf(' %f times, snr = %f , zf error num = %f , mmse error num = %f .\n',looptime,snr,errornum_zf,errornum_mmse);
%            fprintf(' %f times, snr = %f , zf error rate = %f , mmse error rate = %f .\n',looptime,snr,ser_zf,ser_mmse);
        end
            
    end  
    %%
%     ser_ls = gather(ser_ls);

    %%
    if amp == amp_begin
        fsnr_beforecorrect = fopen(save_path+"/snr_beforecorrect.txt",'w');
        fser_beforecorrect = fopen(save_path+"/ser_beforecorrect.txt",'w');
        fsnr_aftercorrect = fopen(save_path+"/snr_aftercorrect.txt",'w');
        fser_aftercorrect = fopen(save_path+"/ser_aftercorrect.txt",'w');
        
        fprintf(fsnr_beforecorrect,' 4pam ,add zero ,pilot length  = %d , zero length  = %d , data length  = %d , bw = %.6g \r\n',pilot_length,zero_length,data_length,origin_rate);
        fprintf(fser_beforecorrect,' 4pam ,add zero ,pilot length  = %d , zero length  = %d , data length  = %d , bw = %.6g \r\n',pilot_length,zero_length,data_length,origin_rate);
        fprintf(fsnr_aftercorrect,' 4pam ,add zero ,pilot length  = %d , zero length  = %d , data length  = %d , bw = %.6g \r\n',pilot_length,zero_length,data_length,origin_rate);
        fprintf(fser_aftercorrect,' 4pam ,add zero ,pilot length  = %d , zero length  = %d , data length  = %d , bw = %.6g \r\n',pilot_length,zero_length,data_length,origin_rate);
    else
        fsnr_beforecorrect = fopen(save_path+"/snr_beforecorrect.txt",'a');
        fser_beforecorrect = fopen(save_path+"/ser_beforecorrect.txt",'a');
        fsnr_aftercorrect = fopen(save_path+"/snr_aftercorrect.txt",'a');
        fser_aftercorrect = fopen(save_path+"/ser_aftercorrect.txt",'a');
    end
    fprintf(fsnr_beforecorrect,'%.8f \r\n',snr_ls_beforecorrect);
    fprintf(fser_beforecorrect,'%.6g \r\n',ser_ls_beforecorrect);
    fprintf(fsnr_aftercorrect,'%.8f \r\n',snr_ls);
    fprintf(fser_aftercorrect,'%.6g \r\n',ser_ls);
    fclose(fsnr_aftercorrect);
    fclose(fser_aftercorrect);

end


