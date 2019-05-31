% example code for large-scale simulation
dendvar = 0.1;
thr_dend = 0.05;
dendinh_width = pi/8;
soma_inh_strs = 1.4:0.2:1.8;
dend_inh_strs = 21:3:27;
intensities = [0.1:0.1:0.8];
contrasts = 0:0.5:2;
num_trial = 50;
for k_soma_inh_str = 1:length(soma_inh_strs)
    for k_dend_inh_str = 1:length(dend_inh_strs)
        for k_intensity = 1:length(intensities)
            for k_contrast = 1:length(contrasts)
                Data{k_soma_inh_str}{k_dend_inh_str}{k_intensity}{k_contrast} = zeros(num_trial,601); % 601 is 'n * 2'(output1, output2) and 401 for soas
                for k_trial = 1:num_trial
                    fprintf('%d (/3) - %d (/3) - %d (/8) - %d (/5) - %d (/50)\n', k_soma_inh_str, k_dend_inh_str, k_intensity, k_contrast, k_trial);
                    [stim, x0, out, soas] = SWMdend1_wofig([100 100],[0 200 400],[pi/12 0 1],[pi/12 0 15],...
                        soma_inh_strs(k_soma_inh_str),[dendinh_width 0 dend_inh_strs(k_dend_inh_str) thr_dend],...
                        [0.1 5 dendvar],[pi/12 intensities(k_intensity) intensities(k_intensity)*contrasts(k_contrast) 0.1],0.05,[]);
                    Data{k_soma_inh_str}{k_dend_inh_str}{k_intensity}{k_contrast}(k_trial,:) = [out(1,:) out(2,:) soas'];
                end
                save Data_sim1 Data
            end
        end
    end
end
