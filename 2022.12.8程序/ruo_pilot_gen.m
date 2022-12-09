% idinput() can also be used to generate m-sequences
function[pilot] = ruo_pilot_gen(pri_poly)
    n = length(pri_poly); % 移位寄存器长�?
    N = 2^n-1; % 伪随机码的周�?
    register = [zeros(1,n-1) 1];
    pilot = zeros(1,N);
    for i = 1:N
        newregister = mod(sum(pri_poly.*register),2);
        register(2:n) = register(1:(n-1));
        register(1) = newregister;
        pilot(i) = register(n);
    end
end

% function[pilot] = pilot_gen(pri_poly)
%     n = length(pri_poly); % 移位寄存器长�?
%     N = 2^n-1; % 伪随机码的周�?
%     register = [zeros(1,n-1) 1];
%     newregister = zeros(1,n);
%     pilot = zeros(1,N);
%     for i = 1:N
%         newregister(1) = mod(sum(pri_poly.*register),2);
%         for j = 2:n
%             newregister(j) = register(j-1);
%         end
%         register = newregister;
%         pilot(i) = register(n);
%     end
% end

% function[pilot] = pilot_gen(pri_poly)
%     n = length(pri_poly); % 移位寄存器长�?
%     N = 2^n-1; % 伪随机码的周�?
%     register = [1 zeros(1,n-2) 1];
%     newregister = zeros(1,n);
%     pilot = zeros(1,N);
%     for i = 1:N
%         pilot(i) = register(n);
%         newregister(1) = mod(sum(pri_poly.*register),2);
%         for j = 2:n
%             newregister(j) = register(j-1);
%         end
%         register = newregister;
%     end
% end

