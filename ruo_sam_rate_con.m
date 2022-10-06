%%%%%%%%%%%%%%%%%  Sample rate conversion %%%%%%%%%%%%%%%

function signal_final_sampled = ruo_sam_rate_con(signal_origin,filter,upf,dof) 
% function signal_final_sampled = sam_rate_con(signal_origin,h,signal_bw,orisam_rate,sample_rate,varargin) 
%signal_bw:signal bandwidth; orisam_rate:signal's origin sampling rate; sample_rate:sampling rate after conversion

%     p = inputParser;
%     p.addParameter('model','transmit',@(x)any(validatestring(x,{'transmit','receive'})));
%     p.addParameter('receive_bandwidth','multiple',@(x)any(validatestring(x,{'signal_bandwidth','multiple'})));
%     p.parse(varargin{:});
%    
% %     default_bw = signal_bw;
%     default_bw = 12.5e6;
% %     multiple_bw = orisam_rate/2.5;
%     multiple_bw = 12.5e6;
% 
%     filter_ord = 1000;     % filter order
%     rp = 0.00057565;
%     rst = 1e-4;
% 
%     llcm = lcm(orisam_rate,sample_rate);
%     upf = llcm/orisam_rate;  % upsampling parameters
%     dof = llcm/sample_rate;   % downsampling parameters
%     ups_rate = orisam_rate*upf;
% 
%     if signal_bw >= ups_rate/2
%         ups_rate = ups_rate*3;
%         upf = ups_rate/orisam_rate;
%         dof = ups_rate/sample_rate;
%     end    

    signal_upsample = upsample(signal_origin,upf);  % upsample signal
    signal_upsample = signal_upsample*upf;

%     switch(p.Results.model)
%         case 'transmit'
%             h = firceqrip(filter_ord,default_bw/(ups_rate/2),[rp rst],'passedge'); %Impulse function
%             h = h./norm(h,2)*sqrt(default_bw/ups_rate*2);
%         case 'receive'
%             switch(p.Results.receive_bandwidth)
%                 case 'signal_bandwidth'
%                     h = firceqrip(filter_ord,default_bw/(ups_rate/2),[rp rst],'passedge'); %Impulse function
%                     h = h./norm(h,2)*sqrt(default_bw/ups_rate*2);
%                 case 'multiple'
%                     h = firceqrip(filter_ord,multiple_bw/(ups_rate/2),[rp rst],'passedge'); %Impulse function
%                     h = h./norm(h,2)*sqrt(multiple_bw/ups_rate*2);
%             end
%     end

    signal_passfilter = conv(signal_upsample,filter);                        % Filter the upsampled signal
    signal_passfilter = signal_passfilter((length(filter)+1)/2:length(signal_passfilter)-(length(filter)-1)/2);
%     signal_passfilter = signal_passfilter./norm(signal_passfilter,2)*sqrt(length(signal_passfilter))*sqrt(bandpower(signal_origin));

    if fix(length(signal_passfilter)/2)+500*dof <= length(signal_passfilter)
        signal_fil_judge = signal_passfilter(1,fix(length(signal_passfilter)/2):fix(length(signal_passfilter)/2)+500*dof);
    else
        signal_fil_judge = signal_passfilter(1,fix(length(signal_passfilter)/2):end);
    end
    energy = zeros(1,dof);
    for downsample_phase = 0:dof-1
        sampled_signal_judge = downsample(signal_fil_judge,dof,downsample_phase);
        energy(downsample_phase+1) = norm(sampled_signal_judge,2);
    end
    [~,judge_phase] = max(energy);
    temporary_phase = fix(length(signal_passfilter)/2)+judge_phase-1;
    real_phase = mod(temporary_phase,dof);
    if real_phase == 0
        real_phase = dof;
    end                                                  % Find the downsample phase with max energy     

    signal_final_sampled = downsample(signal_passfilter,dof,real_phase-1);
%     signal_final_sampled = signal_final_sampled./norm(signal_final_sampled,2)*sqrt(length(signal_final_sampled))*sqrt(bandpower(signal_origin));
    
end


% function sampled_signal = sam_rate_con(signal,signal_bw,orisam_rate,sample_rate) 
% %signal_bw:signal bandwidth; orisam_rate:signal's origin sampling rate; sample_rate:sampling rate after conversion
% 
%     default_bw = 25e6;
%     multiple_bw = 75e6;
% 
%     fil_ord = 1000;     % filter order
%     rp = 0.00057565;
%     rst = 1e-4;
% 
%     llcm = lcm(orisam_rate,sample_rate);
%     upf = llcm/orisam_rate;  % upsampling parameters
%     dof = llcm/sample_rate;   % downsampling parameters
%     ups_rate = orisam_rate*upf;
% 
%     if signal_bw >= ups_rate/2
%         ups_rate = ups_rate*3;
%         upf = ups_rate/orisam_rate;
%         dof = ups_rate/sample_rate;
%     end    
% 
%     sys_us = upsample(signal,upf);  
%     h = firceqrip(fil_ord,default_bw/(ups_rate/2),[rp rst],'passedge');   % Impulse function
%     h = h./norm(h,2);
%     sys_fil = conv(sys_us,h); 
% 
%     sys_fil_judge = sys_fil(1,fix(length(sys_fil)/2):fix(length(sys_fil)/2)+1000*upf);
%     energy = zeros(1,dof);
%     for downsample_phase = 0:dof-1
%         sampled_signal_judge = downsample(sys_fil_judge,dof,downsample_phase);
%         energy(downsample_phase+1) = norm(sampled_signal_judge,2);
%     end
%     [max_energy,judge_phase] = max(energy);
%     temporary_phase = fix(length(sys_fil)/2)+judge_phase-1;
%     real_phase = mod(temporary_phase,dof);
%     if real_phase =  = 0
%         real_phase = dof;
%     end
% 
%     sampled_signal = downsample(sys_fil,dof,real_phase-1);
% end
% 
