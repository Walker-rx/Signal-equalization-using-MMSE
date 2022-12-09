function [fin_syn_point_recorrect] = ruo_signal_syn_recorrect(origin_signal,current_signal,fin_syn_point_ori) 
    times = 1;
    maxlag = 500*times;
    dev = 20;

    signal_ori = origin_signal(1:times*1000); 
    if (fin_syn_point_ori-dev) > 0
        signal_resam = current_signal(fin_syn_point_ori-dev:end);   
    else
        signal_resam = current_signal(1:end);  
    end
             
    if times*3000 <= length(signal_resam)
        signal_rec = signal_resam(1:times*3000);
    else
        signal_rec = signal_resam(1:end);
    end

    [r0,~] = xcorr(signal_rec,signal_ori,maxlag);    
    r1 = abs(r0(maxlag:end));
    [~,index0] = max(r1);                            % Find the point with the largest correlation coefficient
    if index0  ==  1
        fin_syn_point_recorrect = 1;
    else
        fin_syn_point_recorrect = index0-1+fin_syn_point_ori-(dev + 1);                        
    end
    
end


