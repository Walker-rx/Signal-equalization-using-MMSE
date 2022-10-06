clear,close all

% channel_choice = [1];
channel_choice = [4];
quanti = 32000;
dir_up = "./data_set_final/";
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
% for current=780:40:820
%     bias_name = current
    bias_name = 780
%     fprintf(dp821A,strcat(':APPL CH2,8,',num2str(current/1000)));
    pause(1.5);
%     fprintf(dp821A,':MEAS:CURR? CH2');
%     aa = fscanf(dp821A);
%     bias = str2num(aa)
    pause(0.5);
    %%
%     uniform_gen

    %%
%     pam2_gen

    %%
%     pam8_gen

    %%
    pam4_gen

    %%
    % pam16_gen

    %%
%     mseq_gen
% end
% fprintf(dp821A,':OUTP CH2,OFF');
% fprintf(dp821A,':OUTP CH1,OFF');
% fclose(dp821A);



