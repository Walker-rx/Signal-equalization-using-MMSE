% fclose(dp811A);
dp811A = visa('ni','USB0::0x1AB1::0x0E11::DP8D163650047::INSTR');

% dp821A = visa('ni','ASRL3:INSTR');
fopen(dp811A);

current = 300;
fprintf(dp811A,strcat(':APPL CH1,6,',num2str(current/1000)));
fprintf(dp811A,':OUTP CH1,ON');

pause(1.5);
fprintf(dp811A,':MEAS:CURR? CH1');
bias = str2num(fscanf(dp811A))

fprintf(dp811A,':OUTP CH1,OFF');
fclose(dp811A);