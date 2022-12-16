clear;
close all;

channel_choice = 1; 
channel_choice_inf = 3;
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

M = 8;
data_length_initial = 10000;
zero_length = 3000;
zero_length_forsyn = 200;
ls_order = 50;
num_of_windows = 100;

times = 6;
origin_rate_tmp = 10e6;
f_rate = 160e6;
d_rate_tmp = origin_rate_tmp*times;

filter_order = 1000;     % filter order used in function sam_rate_con
rp = 0.00057565;      
rst = 1e-4;       % filter parameter used in function sam_rate_con

%% Generate pilot
% pilot = pilot_gen([0 0 0 1 1 1 1]);
% pilot = pilot_gen([1 1 1 0 0 0 1 1 1]);
pilot = ruo_pilot_gen([1 1 1 1 1 1 0 1 0 0 1]);    %  2^11 = 2048
% pilot = ruo_pilot_gen([1 1 1 1 1 1 0 1 0 1 0 1 1]);   %  2^13 = 8192
pilot_ps = bandpower(pilot);
pilot_length = length(pilot);
pilot_bpsk = pilot*2-1;

%% Generate filter
[origin_rate , d_rate , upf_transmit , dof_transmit , filter_transmit , upf_receive , dof_receive , filter_receive] = ruo_filter_gen(origin_rate_tmp , f_rate , times , filter_order , rp , rst);
rate_change = 100*abs(origin_rate-origin_rate_tmp)/origin_rate_tmp;
fprintf("rate's change = %.8f%% \n",rate_change);
ups_time = upf_transmit/dof_transmit;
data_length = round(data_length_initial*8/ups_time);
%%
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

%%
t = datetime('now');
save_path = "snr_ser/direct/not_replace/"+t.Month+"."+t.Day+"/"+origin_rate_tmp/1e6+"M/"+M+"pam";
if(~exist(save_path,'dir'))
    mkdir(char(save_path));
end

amp_begin = -8;
amp_end = 90;
amp_inf = 40;
fprintf('add zero,ls order=%d,pilot length=%d .\n',ls_order,pilot_length);

for amp = amp_begin:amp_end
    tic
    fprintf("amp = %d \n",amp);
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

    while(errornum_ls_aftercorrect <= 100 || looptime < 300)
%     while(looptime < 10)
        
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
      
        ruo_send;

        %% Calculating the SER without correction
        signal_received_real_send1 = ruo_sam_rate_con(signal_pass_channel,filter_receive,upf_receive,dof_receive);       
        [fin_syn_point_real_send1_tmp , coar_syn_point_real_send1_tmp] = ruo_signal_syn(origin_rate,d_rate,signal_ori,signal_received_real_send1,num_of_windows);      
        if (length(signal_received_real_send1) - length(signal_received_real_send1(fin_syn_point_real_send1_tmp:end))) > (pilot_length + zero_length_forsyn)*times
            coar_syn_point_real_send1 = coar_syn_point_real_send1_tmp;
            fin_syn_point_real_send1 = fin_syn_point_real_send1_tmp;
        else
            coar_syn_point_real_send1 = coar_syn_point_real_send1_tmp + pilot_length + zero_length_forsyn;
            fin_syn_point_real_send1 = fin_syn_point_real_send1_tmp + (pilot_length + zero_length_forsyn)*times;
        end              
        [length_loop_real_send1 , ps_loop_real_send1 , pn_loop_real_send1 , ...
            errornum_ls_loop_real_send1 , error_location_loop_real_send1 , data_demod_ls_real_send1] ...
            = ruo_calculate_ser( data , signal_ori , signal_received_real_send1 , pilot_length , zero_length , data_length ...
                                 , fin_syn_point_real_send1 , times , ls_order);
                                           
        total_length_beforecorrect = total_length_beforecorrect + length_loop_real_send1 ;
        ps_beforecorrect = ps_beforecorrect + ps_loop_real_send1 ;
        pn_beforecorrect = pn_beforecorrect + pn_loop_real_send1 ;
        errornum_ls_beforecorrect = errornum_ls_beforecorrect + errornum_ls_loop_real_send1 ;
        
        snr_ls_beforecorrect = 10*log10(ps_beforecorrect/pn_beforecorrect);
        ser_ls_beforecorrect = errornum_ls_beforecorrect/total_length_beforecorrect;
        
        fprintf(' amp = %d , %f times, before error = %d  .\n',amp,looptime , errornum_ls_loop_real_send1 );
        if mod(looptime,5) == 0
           fprintf('\n');
           fprintf(' amp = %d , %f times, snr = %d , total num = %d , total ls error num = %d , ls error rate = %.6g .\n',amp,looptime,snr_ls_beforecorrect,total_length_beforecorrect,errornum_ls_beforecorrect,ser_ls_beforecorrect);
%            pause(1)
        end
            
    end  

    %%

    %% Save data               
    if amp == amp_begin
        fsnr_beforecorrect = fopen(save_path+"/snr_beforecorrect.txt",'w');
        fser_beforecorrect = fopen(save_path+"/ser_beforecorrect.txt",'w');
        
        fprintf(fsnr_beforecorrect,' %dpam ,add zero ,pilot length  = %d , zero length  = %d , data length  = %d , origin_rate = %.6g \r\n',M,pilot_length,zero_length,data_length,origin_rate);
        fprintf(fser_beforecorrect,' %dpam ,add zero ,pilot length  = %d , zero length  = %d , data length  = %d , origin_rate = %.6g \r\n',M,pilot_length,zero_length,data_length,origin_rate);
    else
        fsnr_beforecorrect = fopen(save_path+"/snr_beforecorrect.txt",'a');
        fser_beforecorrect = fopen(save_path+"/ser_beforecorrect.txt",'a');
    end
    fprintf(fsnr_beforecorrect,'%.8f \r\n',snr_ls_beforecorrect);
    fprintf(fser_beforecorrect,'%.6g \r\n',ser_ls_beforecorrect);
   
    fclose(fsnr_beforecorrect);
    fclose(fser_beforecorrect);
    toc
end


