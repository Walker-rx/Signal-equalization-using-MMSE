clear;
close all;

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
zero_length_forsyn = 200;
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

% filter_ord = 1000;     % filter order used in function sam_rate_con
filter_ord = 200;     % filter order used in function sam_rate_con
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

t = datetime('now');
% save_path = "snr_ser/direct/"+t.Month+"."+t.Day+"/correct_transmit";
save_path = "snr_ser/direct/12.5/correct_transmit";
if(~exist(save_path,'dir'))
    mkdir(char(save_path));
end

amp_begin = -6;
amp_end = 62;
amp_inf = 33;
fprintf('add zero,ls order=%d,pilot length=%d .\n',ls_order,pilot_length);

for amp = 35:amp_end
    looptime = 0;
    ps_aftercorrect = 0;
    pn_aftercorrect = 0;
    ps_beforecorrect = 0;
    pn_beforecorrect = 0;
    errornum_zf = 0;
    errornum_mmse = 0;
    errornum_ls_aftercorrect = 0;
    errornum_ls_beforecorrect = 0;
    total_length_beforecorrect = 0;
    total_length_aftercorrect = 0;
    txerror_num = 0;
    replace_length = 8*times;
    replace_valid_num = 0;
    replace_correct_num = 0;
    fprintf('amp = %d .\n', amp);
    
    while(errornum_ls_aftercorrect <= 30 || looptime < 1000)
%     while(errornum_zf <= 100 || errornum_mmse <= 100 || looptime < 50)
        
        looptime = looptime+1;
         %% Signal send
        data = randi([0,M-1],[1,data_length]);
        data_mpam = real(pammod(data,M));
        data_mpam_ps = bandpower(data_mpam);  
        
        unmod = [pilot,zeros(1,zero_length),data];
        signal_ori = [pilot_bpsk,zeros(1,zero_length),data_mpam];
        pilot_bpsk_forsyn = [pilot_bpsk,zeros(1,zero_length_forsyn)];
        
        pilot_upsample_forsyn = ruo_sam_rate_con(pilot_bpsk_forsyn,filter_transmit,upf_transmit,dof_transmit);
        signal_upsample = ruo_sam_rate_con(signal_ori,filter_transmit,upf_transmit,dof_transmit);

        pilot_send_forsyn = pilot_upsample_forsyn./norm(pilot_upsample_forsyn,2)*sqrt(length(pilot_upsample_forsyn))*100*1.1^(amp_inf);
        signal_send_tmp = signal_upsample./norm(signal_upsample,2)*sqrt(length(signal_upsample))*100*1.1^(amp);
        signal_send_tmp_inf = signal_upsample./norm(signal_upsample,2)*sqrt(length(signal_upsample))*100*1.1^(amp_inf);
       
        signal_rand = rand([1,512]);
        signal_send = [signal_rand pilot_send_forsyn signal_send_tmp];
        signal_send_inf = [signal_rand pilot_send_forsyn signal_send_tmp_inf];
      
        ruo_pam4_send;
        ruo_pam4_send_correct;
        
        %% Locating the transmission error point
        signal_received_inf = ruo_sam_rate_con(signal_pass_channel_inf,filter_receive,upf_receive,dof_receive);
        
        [fin_syn_point_inf_tmp,coar_syn_point_inf_tmp] = ruo_signal_syn(origin_rate,d_rate,signal_ori,signal_received_inf,num_of_windows);
        
        coar_syn_point_inf = coar_syn_point_inf_tmp + pilot_length + zero_length_forsyn;
        fin_syn_point_inf = fin_syn_point_inf_tmp + (pilot_length + zero_length_forsyn)*times;
        
        signal_demod_ls_inf = ruo_signal_equal_ls(signal_ori,signal_received_inf,times,fin_syn_point_inf,pilot_length,zero_length,ls_order);
  
        data_demod_ls_inf = signal_demod_ls_inf(pilot_length+zero_length+1:end);
        
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
        signal_received_channel1_send1 = ruo_sam_rate_con(signal_pass_channel,filter_receive,upf_receive,dof_receive);       
        [fin_syn_point_channel1_send1_tmp , coar_syn_point_channel1_send1_tmp] = ruo_signal_syn(origin_rate,d_rate,signal_ori,signal_received_channel1_send1,num_of_windows);      
        if (length(signal_received_channel1_send1) - length(signal_received_channel1_send1(fin_syn_point_channel1_send1_tmp:end))) > (pilot_length + zero_length_forsyn)*times
            coar_syn_point_channel1_send1 = coar_syn_point_channel1_send1_tmp;
            fin_syn_point_channel1_send1 = fin_syn_point_channel1_send1_tmp;
        else
            coar_syn_point_channel1_send1 = coar_syn_point_channel1_send1_tmp + pilot_length + zero_length_forsyn;
            fin_syn_point_channel1_send1 = fin_syn_point_channel1_send1_tmp + (pilot_length + zero_length_forsyn)*times;
        end              
        [length_loop_channel1_send1 , ps_loop_channel1_send1 , pn_loop_channel1_send1 , ...
            errornum_ls_loop_channel1_send1 , error_location_loop_channel1_send1 , data_demod_ls_channel1_send1] ...
            = ruo_calculate_ser( data , signal_ori , signal_received_channel1_send1 , pilot_length , zero_length , data_length ...
                                 , fin_syn_point_channel1_send1 , times , ls_order);
                                   
        signal_received_channel1_send2 = ruo_sam_rate_con(signal_pass_channel_correct,filter_receive,upf_receive,dof_receive);
        [fin_syn_point_channel1_send2_tmp , coar_syn_point_channel1_send2_tmp] = ruo_signal_syn(origin_rate,d_rate,signal_ori,signal_received_channel1_send2,num_of_windows);
        if (length(signal_received_channel1_send2) - length(signal_received_channel1_send2(fin_syn_point_channel1_send2_tmp:end))) > (pilot_length + zero_length_forsyn)*times
            coar_syn_point_channel1_send2 = coar_syn_point_channel1_send2_tmp;
            fin_syn_point_channel1_send2 = fin_syn_point_channel1_send2_tmp;
        else
            coar_syn_point_channel1_send2 = coar_syn_point_channel1_send2_tmp + pilot_length + zero_length_forsyn;
            fin_syn_point_channel1_send2 = fin_syn_point_channel1_send2_tmp + (pilot_length + zero_length_forsyn)*times;
        end        
        [length_loop_channel1_send2 , ps_loop_channel1_send2 , pn_loop_channel1_send2 , ...
            errornum_ls_loop_channel1_send2 , error_location_loop_channel1_send2 , data_demod_ls_channel1_send2] ...
            = ruo_calculate_ser( data , signal_ori , signal_received_channel1_send2 , pilot_length , zero_length , data_length ...
                                 , fin_syn_point_channel1_send2 , times , ls_order);
        
        total_length_beforecorrect = total_length_beforecorrect + length_loop_channel1_send1 + length_loop_channel1_send2;
        ps_beforecorrect = ps_beforecorrect + ps_loop_channel1_send1 + ps_loop_channel1_send2;
        pn_beforecorrect = pn_beforecorrect + pn_loop_channel1_send1 + pn_loop_channel1_send2;
        errornum_ls_beforecorrect = errornum_ls_beforecorrect + errornum_ls_loop_channel1_send1 + errornum_ls_loop_channel1_send2;
        
        snr_ls_beforecorrect = 10*log10(ps_beforecorrect/pn_beforecorrect);
        ser_ls_beforecorrect = errornum_ls_beforecorrect/total_length_beforecorrect;
        
        %% Replacing the transmission error point
        signal_received_tmp = signal_received_channel1_send1;
        signal_received_aftercorrect = signal_received_channel1_send1;
        signal_ori_resyn = signal_received_aftercorrect(fin_syn_point_channel1_send1_tmp:end);
        fin_syn_point_forcorrect_tmp = ruo_signal_syn_recorrect(signal_ori_resyn,signal_received_channel1_send2,fin_syn_point_channel1_send2_tmp);
        if fin_syn_point_channel1_send1 == fin_syn_point_channel1_send1_tmp && fin_syn_point_channel1_send2 == fin_syn_point_channel1_send2_tmp
            fin_syn_point_forcorrect = fin_syn_point_forcorrect_tmp;
        else         
            fin_syn_point_forcorrect = fin_syn_point_forcorrect_tmp + (pilot_length + zero_length_forsyn)*times;
        end
        if txerror_num_loop > 0
            replace_loc = error_location_inf*times + fin_syn_point_channel1_send1 - 1;
            replace_loc_correct = error_location_inf*times + fin_syn_point_forcorrect - 1;         
            for i = 1:txerror_num_loop
                if (replace_loc(i)-replace_length) > 0
                    if replace_loc(i)+replace_length > length(signal_received_aftercorrect) || replace_loc_correct(i)+replace_length > length(signal_received_channel1_send2)
                        signal_received_aftercorrect(length(signal_received_aftercorrect)-2*replace_length : end) = signal_received_channel1_send2(length(signal_received_channel1_send2)-2*replace_length : end);
                    else
                        signal_received_aftercorrect(replace_loc(i)-replace_length : replace_loc(i)+replace_length) = signal_received_channel1_send2(replace_loc_correct(i)-replace_length : replace_loc_correct(i)+replace_length);
                    end                   
                else
                    signal_received_aftercorrect(1 : replace_loc(i)+replace_length) = signal_received_channel1_send2(1 : replace_loc_correct(i)+replace_length);
                end
            end
        end
        %% Calculating the SER after correction
        [length_loop_aftercorrect , ps_loop_aftercorrect , pn_loop_aftercorrect , ...
            errornum_ls_loop_aftercorrect , error_location_aftercorrect , data_demod_ls_aftercorrect] ...
            = ruo_calculate_ser( data , signal_ori , signal_received_aftercorrect , pilot_length , zero_length , data_length ...
                                 , fin_syn_point_channel1_send1 , times , ls_order);
        
        total_length_aftercorrect = total_length_aftercorrect + length_loop_aftercorrect;
        errornum_ls_aftercorrect = errornum_ls_aftercorrect + errornum_ls_loop_aftercorrect;
        ps_aftercorrect = ps_aftercorrect + ps_loop_aftercorrect;
        pn_aftercorrect = pn_aftercorrect + pn_loop_aftercorrect;
        
        ser_ls_aftercorrect = errornum_ls_aftercorrect/total_length_aftercorrect;
        snr_ls_aftercorrect = 10*log10(ps_aftercorrect/pn_aftercorrect);
        
        replace_correct_num_tmp = errornum_ls_loop_channel1_send1 - errornum_ls_loop_aftercorrect;
        replace_correct_num = replace_correct_num + replace_correct_num_tmp;
        if errornum_ls_loop_channel1_send1 > errornum_ls_loop_aftercorrect
            replace_valid_num = replace_valid_num + 1;
        end
%         snr_ls_beforecorrect = 0;
%         ser_ls_beforecorrect = 0; 
        if mod(looptime,5) == 0
           fprintf(' amp = %d , %d times, average tx error num = %d .\n',amp,looptime,txerror_num_average);
           fprintf(' %f times, before, snr = %d , total ls error num = %d , ls error rate = %.6g .\n',looptime,snr_ls_beforecorrect,errornum_ls_beforecorrect/2,ser_ls_beforecorrect);
           fprintf(' %f times, after, snr = %d , total ls error num = %d , ls error rate = %.6g .\n',looptime,snr_ls_aftercorrect,errornum_ls_aftercorrect,ser_ls_aftercorrect);
           fprintf(' %f times, before error = %d , replace num = %d , after error = %d .\n',looptime , errornum_ls_loop_channel1_send1 , txerror_num_loop ,  errornum_ls_loop_aftercorrect);
%            disp(["error location = ",error_location]);
%            disp(["correct = ",data(error_location)]);
%            disp(["false = ",data_demod_ls(error_location)]);
%            fprintf(' %f times, snr = %f , zf error num = %f , mmse error num = %f .\n',looptime,snr,errornum_zf,errornum_mmse);
%            fprintf(' %f times, snr = %f , zf error rate = %f , mmse error rate = %f .\n',looptime,snr,ser_zf,ser_mmse);
        end
            
    end  
    %%
%     ser_ls = gather(ser_ls);

    %% Save data               
    if amp == amp_begin
        fsnr_beforecorrect = fopen(save_path+"/snr_beforecorrect.txt",'w');
        fser_beforecorrect = fopen(save_path+"/ser_beforecorrect.txt",'w');
        fsnr_aftercorrect = fopen(save_path+"/snr_aftercorrect.txt",'w');
        fser_aftercorrect = fopen(save_path+"/ser_aftercorrect.txt",'w');
        freplace_invalid = fopen(save_path+"/replace_invalid.txt",'w');
        
        fprintf(fsnr_beforecorrect,' 4pam ,add zero ,pilot length  = %d , zero length  = %d , data length  = %d , bw = %.6g \r\n',pilot_length,zero_length,data_length,origin_rate);
        fprintf(fser_beforecorrect,' 4pam ,add zero ,pilot length  = %d , zero length  = %d , data length  = %d , bw = %.6g \r\n',pilot_length,zero_length,data_length,origin_rate);
        fprintf(fsnr_aftercorrect,' 4pam ,add zero ,pilot length  = %d , zero length  = %d , data length  = %d , bw = %.6g \r\n',pilot_length,zero_length,data_length,origin_rate);
        fprintf(fser_aftercorrect,' 4pam ,add zero ,pilot length  = %d , zero length  = %d , data length  = %d , bw = %.6g \r\n',pilot_length,zero_length,data_length,origin_rate);
    else
        fsnr_beforecorrect = fopen(save_path+"/snr_beforecorrect.txt",'a');
        fser_beforecorrect = fopen(save_path+"/ser_beforecorrect.txt",'a');
        fsnr_aftercorrect = fopen(save_path+"/snr_aftercorrect.txt",'a');
        fser_aftercorrect = fopen(save_path+"/ser_aftercorrect.txt",'a');
        freplace_invalid = fopen(save_path+"/replace_invalid.txt",'a');
    end
    fprintf(fsnr_beforecorrect,'%.8f \r\n',snr_ls_beforecorrect);
    fprintf(fser_beforecorrect,'%.6g \r\n',ser_ls_beforecorrect);
    fprintf(fsnr_aftercorrect,'%.8f \r\n',snr_ls_aftercorrect);
    fprintf(fser_aftercorrect,'%.6g \r\n',ser_ls_aftercorrect);
    fprintf(freplace_invalid,'amp = %.d , looptime = %d , replace valid num = %d , replace correct num = %d \r\n',amp,looptime,replace_valid_num,replace_correct_num);

    fclose(fsnr_beforecorrect);
    fclose(fser_beforecorrect);
    fclose(fsnr_aftercorrect);
    fclose(fser_aftercorrect);
    fclose(freplace_invalid);
end


