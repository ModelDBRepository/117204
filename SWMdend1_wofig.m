function [stimulus, x0, outputs, soas] = SWMdend1_wofig(p_mn, timings, p_ff, p_re, q_ris, p_rid, p_dend, p_stim, initrand, isavename)

% Spatial Working Memory Circuit Incorporating Nonlinear Input Integration in the Dendritic Branches
%
%   p_mn : [m n] (number of dendritic branches per neuron / number of neurons) % e.g., [50 50]
%   timings : [t_stim_start, t_stim_end, t_final] % e.g., [0 100 200]
%   p_ff : [Jsigma_ff, Jbase_ff, Jamp_ff] % e.g., [pi/12 0 1]
%   p_re : [Jsigma_re, Jbase_re, Jamp_re] % e.g., [pi/12 0 15]
%   q_ris : str_ris % e.g., 2
%   p_rid : [Jsigma_rid, Jbase_rid, Jamp_rid, thr_rid] % e.g., [pi/2 2 0 0], [pi/4 0 16 0]
%   p_dend : [dend_threshold * m (number of branches), dend_slope, dendrand] % e.g., [0 1 0]
%           #dendrand : randomness on the parameters regarding the dendritic nonlinearity, i.e., (standard variation / mean) of threshold, slope, and saturation
%   p_stim : [Jsigma_stimulus, Jbase_stimulus, Jamp_stimulus, noiselevel] % e.g., [pi/12 0.5 0 0.2] (BLA), or [pi/12 0.5 0.5 0.2] (WTA)
%   initrand : randomness on the initial activities of the neurons will be set to "initrand * rand(n,1)" % e.g., 0.1 or 0 or 0.05
%
%   stimulus : actually presented stimulus
%   x0 : actual initial activites of the neurons
%   outputs: [(1st row) neural activities at the end of the stimulus presentation; (2nd row) neural activities at the end of the simulation]
%   soas : time evolution of the "Sum Of Activity"
%   
%	(e.g)
% dendvar 0.1 + thrdend 0.05 : pi/8 dendritic
% [stim, x0, out, soas] = SWMdend1_wofig([100 100],[0 200 400],[pi/12 0 1],[pi/12 0 15],1.6,[pi/8 0 24 0.05],[0.1 5 0.1],[pi/12 0.3 0.1 0.1],0.05,[]);
%
% Kenji Morita
% Mar 2008 @ Tokyo


%%%%% the input dimension and the number of the principal neurons %%%%%
global m n
m = p_mn(1); % input dimension
n = p_mn(2); % number of the principal neurons


%%%%% protocol %%%%%
global t_stim_start t_stim_end t_final
t_stim_start = timings(1);
t_stim_end = timings(2);
t_final = timings(3);


%%%%% connection strengths %%%%%
global ff re str_ris rid thr_rid

% phase (selectivity) of the input layer ('pre') and the output layer ('post')
phase_pre = 2*pi * [0:m-1] / m;
phase_post = 2*pi * [0:n-1] / n;

% feed-forward excitation from the input layer to the output layer
ff = zeros(n,m);
Jsigma_ff = p_ff(1);
Jbase_ff = p_ff(2); % constant baseline (NB: this is not minimum value of the connection)
Jamp_ff = p_ff(3); % amplitude of the exp-cos curve
tmp = (phase_post' * ones(1,m)) - (ones(n,1) * phase_pre); % self - other
phasedifs_ff = min(abs(tmp), 2*pi - abs(tmp)); % NB: be aware of the necessity of this step
ff = (Jamp_ff * exp((cos(phasedifs_ff) - 1) / Jsigma_ff^2) + Jbase_ff) / m; % NB: should be normalized w.r.t. the input dimension

% recurrent excitation between the principal neurons in the output layer
re = zeros(n,n);
Jsigma_re = p_re(1);
Jbase_re = p_re(2); % constant baseline (NB: this is not minimum value of the connection)
Jamp_re = p_re(3); % amplitude of the exp-cos curve
tmp = (phase_post' * ones(1,n)) - (ones(n,1) * phase_post); % self - other
phasedifs_re = min(abs(tmp), 2*pi - abs(tmp)); % NB: be aware of the necessity of this step
re = (Jamp_re * exp((cos(phasedifs_re) - 1) / Jsigma_re^2) + Jbase_re) / n; % NB: should be normalized w.r.t. the number of the neurons

% recurrent somatic inhibition between the principal neurons in the output layer, via soma-targeting GABAergic cell
str_ris = q_ris; % strength of the somatic inhibition

% recurrent dendritic inhibition between the principal neurons in the output layer, via hidden dendrite-targeting GABAergic cells
rid = zeros(n,n);
Jsigma_rid = p_rid(1);
Jbase_rid = p_rid(2); % constant baseline (NB: this is not minimum value of the connection)
Jamp_rid = p_rid(3); % amplitude of the exp-cos curve
tmp = (phase_post' * ones(1,n)) - (ones(n,1) * phase_post); % self - other
phasedifs_rid = min(abs(tmp), 2*pi - abs(tmp)); % NB: be aware of the necessity of this step
rid = (Jamp_rid * exp((cos(phasedifs_rid) - 1) / Jsigma_rid^2) + Jbase_rid) / n; % NB: should be normalized w.r.t. the number of the neurons
thr_rid = p_rid(4); % threshold for dendritic recurrent inhibition
%plot(rid(round(n)/2+1,:));


%%%%% dendritic nonlinearity %%%%%
global dend_threshold dend_slope dend_upperbound % NB: these are n*m matrice
dend_threshold = p_dend(1) * (1/m) * max((1 + p_dend(3) * randn(n,m)), 0);
dend_slope = p_dend(2) * max((1 + p_dend(3) * randn(n,m)), 0);
dend_upperbound = 1/m * max((1 + p_dend(3) * randn(n,m)), 0);


%%%%% input stimulus %%%%%
global stimulus
global stimulus_signal stimulus_noise_level
stimulus = zeros(1,m);
stimcenter = pi; % location (radian) of the stimulus center
Jsigma_stimulus = p_stim(1);
Jbase_stimulus = p_stim(2); % constant baseline (NB: this is not minimum value of the connection)
Jamp_stimulus = p_stim(3); % amplitude of the exp-cos curve
noiselevel = p_stim(4);
tmp = stimcenter - phase_pre; % self - other
phasedifs_stimulus = min(abs(tmp), 2*pi - abs(tmp)); % NB: be aware of the necessity of this step
stimulus_signal = (Jamp_stimulus * exp((cos(phasedifs_stimulus) - 1) / Jsigma_stimulus^2) + Jbase_stimulus);
stimulus_noise_level = noiselevel * Jbase_stimulus;
%stimulus_noise = noiselevel * Jbase_stimulus * randn(1,m);
%stimulus = stimulus_signal + stimulus_noise;
if 0
    F = figure;
    A = axes;
    set(F,'position',[4    406   506   287]);
    hold on;
    tmp = bar([0.5:1:m-0.5],stimulus_signal,1);
    set(tmp,'BarWidth',1,'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0.5 0.5 0.5]);
    %shading interp
    plot([m/2+0.5 m/2+0.5],[0 1.2],'k:');
    axis([0 m 0 1.2]);
    set(A,'Box','on');
    set(A,'FontName','Ariel','FontSize',18);
	set(A,'XTick',[0 m/2+0.5 m],'XTickLabel',[]);
    set(A,'YTick',[0:0.2:1.2],'YTickLabel',[]);
    if ~isempty(isavename)
        print(F,'-depsc',[isavename '_input']);
    end
end


%%%%% main ode loop %%%%%
global input_ff input_re input_ris input_rid dend_upperbound J_dend % need to refer to those values to plot each component of the input
x0 = initrand * rand(n,1);
if 1 % Runge-Kutta
    tstep_calc = 0.1;
    tstep_save = 1;
    % until the end of the stimulus presentation
    [tseq1, xseq1] = runge4('SWMdend1_fun1', 0, t_stim_end, x0, tstep_calc, tstep_save);
    output1 = xseq1(end,:);
    % after the extinction of the stimulus
    [tseq2, xseq2] = runge4('SWMdend1_fun1', t_stim_end, t_final, output1', tstep_calc, tstep_save);
    output2 = xseq2(end,:);
    % integrate
    tseq = [tseq1(1:end-1); tseq2];
    xseq = [xseq1(1:end-1,:); xseq2];
elseif 0 % ode45
    odeopt = odeset('NonNegative', [1:n], 'MaxStep', 1);
    
    % until the end of the stimulus presentation
    [tseq1,xseq1]=ode45(@SWMdend1_fun1,[0,t_stim_end],x0,odeopt);
    output1 = xseq1(end,:);
    if 0 % plot the input components
        
        % specify the potential winner and the potential loser neurons
        [tmpvalue,pot_win_index] = min(abs(phase_post - stimcenter)); pot_win_index = min(pot_win_index);
        [tmpvalue,pot_los_index] = min(abs(phase_post - mod((stimcenter - pi),2*pi))); pot_los_index = min(pot_los_index);
        
        % for the potential winner
        F = figure;
        A = axes;
        set(F,'Position',[516   406   506   287]);
        hold on
        plot(1:m,input_ff(pot_win_index,:),'m');
        plot([1 m],input_re(pot_win_index)/m * ones(1,2),'r');
        plot(1:m,J_dend(pot_win_index,:),'g');
        plot([1 m],dend_threshold * ones(1,2),'c');
        plot([1 m],dend_threshold+dend_upperbound * ones(1,2),'y');
        plot([1 m],[0 0],'k');
        axis([1 m -dend_upperbound/10 dend_upperbound*1.1]);
        set(A,'Box','on');
        set(A,'FontName','Ariel','FontSize',18);
        set(A,'XTick',[1 (1+m)/2 m],'XTickLabel',[]);
        set(A,'YTick',[0],'YTickLabel',[0]);
        if ~isempty(isavename)
            print(F,'-depsc',[isavename '_pot_win_inputs']);
        end
        
        % for the potential loser
        F = figure;
        A = axes;
        set(F,'Position',[516    37   506   287]);
        hold on
        plot(1:m,input_ff(pot_los_index,:),'m');
        plot([1 m],input_re(pot_los_index)/m * ones(1,2),'r');
        plot(1:m,J_dend(pot_los_index,:),'g');
        plot([1 m],dend_threshold * ones(1,2),'c');
        plot([1 m],dend_threshold+dend_upperbound * ones(1,2),'y');
        plot([1 m],[0 0],'k');
        axis([1 m -dend_upperbound/10 dend_upperbound*1.1]);
        set(A,'Box','on');
        set(A,'FontName','Ariel','FontSize',18);
        set(A,'XTick',[1 (1+m)/2 m],'XTickLabel',[]);
        set(A,'YTick',[0],'YTickLabel',[0]);
        if ~isempty(isavename)
            print(F,'-depsc',[isavename '_pot_los_inputs']);
        end

    end
    
    % after the extinction of the stimulus
    [tseq2,xseq2]=ode45(@SWMdend1_fun1,[t_stim_end,t_final],output1',odeopt);
    output2 = xseq2(end,:);
    % integrate
    tmp_tseq = [tseq1; tseq2(2:end)];
    tmp_xseq = [xseq1; xseq2(2:end,:)];
    tseq = [0:t_final]';
    xseq = interp1(tmp_tseq,tmp_xseq,tseq);
end

% output variables
outputs = [output1; output2];
soas = sum(xseq,2);


%%%%% plot the neural activities %%%%%
% stimulus end
if 0
    % calculate the peak location
	tmp_cos = sum(output1 .* cos(2*pi*[0:1/n:(n-1)/n]));
    tmp_sin = sum(output1 .* sin(2*pi*[0:1/n:(n-1)/n]));
	tmp_peak = atan(tmp_sin / tmp_cos) / (2*pi);
	if tmp_cos < 0
        tmp_peak = tmp_peak + 0.5;
	elseif tmp_sin < 0 % (& tmp_cos > 0)
        tmp_peak = tmp_peak + 1;
    end
    peak1 = tmp_peak;
	% plot
    F = figure;
    A = axes;
    set(F,'Position',[516   406   506   287]);
    hold on;
    tmp = bar([0.5:1:n-0.5],output1,1);
    set(tmp,'BarWidth',1,'FaceColor','r','Edgecolor','r');
    plot([n/2+0.5 n/2+0.5],[0 1.2],'k:');
    plot([n*peak1+0.5 n*peak1+0.5],[0 1.2],'k'); % NB: 'n*' is necessary, and also, '+0.5' is necessary to correct the shift of 'bar'
    axis([0 n 0 1.2]); % for 'bar'
    set(A,'Box','on');
    set(A,'FontName','Ariel','FontSize',18);
    set(A,'XTick',[0 n/2+0.5 n],'XTickLabel',[]); % for 'bar'
    set(A,'YTick',[0 1],'YTickLabel',[]);
    %set(A,'PlotBoxAspectRatio',[2 1 1]);
    if ~isempty(isavename)
        print(F,'-depsc',[isavename '_output1']);
    end
end
% final
if 0
    % calculate the peak location
	tmp_cos = sum(output2 .* cos(2*pi*[0:1/n:(n-1)/n]));
    tmp_sin = sum(output2 .* sin(2*pi*[0:1/n:(n-1)/n]));
	tmp_peak = atan(tmp_sin / tmp_cos) / (2*pi);
	if tmp_cos < 0
        tmp_peak = tmp_peak + 0.5;
	elseif tmp_sin < 0 % (& tmp_cos > 0)
        tmp_peak = tmp_peak + 1;
    end
    peak2 = tmp_peak;
    % plot
    F = figure;
    A = axes;
    set(F,'Position',[516    37   506   287]);
    hold on;
    tmp = bar([0.5:1:n-0.5],output2,1);
    set(tmp,'BarWidth',1,'FaceColor','r','Edgecolor','r');
    plot([n/2+0.5 n/2+0.5],[0 1.2],'k:');
    plot([n*peak2+0.5 n*peak2+0.5],[0 1.2],'k'); % NB: 'n*' is necessary, and also, '+0.5' is necessary to correct the shift of 'bar'
    axis([0 n 0 1.2]); % for 'bar'
    set(A,'Box','on');
    set(A,'FontName','Ariel','FontSize',18);
    set(A,'XTick',[0 n/2+0.5 n],'XTickLabel',[]); % for 'bar'
    set(A,'YTick',[0 1],'YTickLabel',[]);
    %set(A,'PlotBoxAspectRatio',[2 1 1]);
    if ~isempty(isavename)
        print(F,'-depsc',[isavename '_output2']);
    end
end
% activity evolution 2D
if 0
    F = figure;
    A = axes;
    set(F,'Position',[4    37   506   287]);
    hold on;
    image_amp = 60; % just for nice view
    image_power = 1/3; % just for nice view
    image_matrix = max(0, fliplr(xseq')'); % just to avoid infinitesimal negative value "-0.0000", which causes problem if imposed by '.^(1/3) etc'.
    O = image(image_amp * image_matrix.^(image_power));
    L = plot([0.5+(n/2+0.5) 0.5+(n/2+0.5)],[0.5 t_final+0.5],'w'); % '0.5+' is because 'image' plot starts at 0.5. '+0.5)' is because of the band width.
    %C = colorbar;
    axis([0.5 n+0.5 0.5 t_final+1.5]); % '0.5' is because 'image' plot starts at 0.5.
    set(A,'Box','on');
    set(A,'FontName','Ariel','FontSize',18);
    set(A,'XTick',[0.5 0.5+(n/2+0.5) n+0.5],'XTickLabel',[]);
    set(A,'YTick',[0.5 200.5 401.5],'YTickLabel',[400 200 0]);
    set(A,'PlotBoxAspectRatio',[1 1.5 1]);
    %set(C,'YTick',[0 image_amp],'YTickLabel',[]);
    if ~isempty(isavename)
        print(F,'-depsc',[isavename '_2D']);
    end
end



%%%%% Runge-Kutta %%%%%
function [tseq, xseq] = runge4(fun_name, t_init, t_end, x_init, tstep_calc, tstep_save)

howoftensave = tstep_save / tstep_calc; % this SHOULD BE an integer
howmanystep = (t_end - t_init) / tstep_calc; % this SHOULD BE an integer

tseq = [t_init:tstep_save:t_end]';
xseq = zeros(length(tseq),length(x_init));
xseq(1,:) = x_init';

t = t_init;
x = x_init;

for k = 1:howmanystep
    tmp1 = tstep_calc * feval(fun_name, t, x);
    tmp2 = tstep_calc * feval(fun_name, t+tstep_calc/2, x+tmp1/2);
    tmp3 = tstep_calc * feval(fun_name, t+tstep_calc/2, x+tmp2/2);
    tmp4 = tstep_calc * feval(fun_name, t+tstep_calc, x+tmp3);
    delta = (tmp1 + tmp4) / 6 + (tmp2 + tmp3) / 3;
    x = x + delta;
    if mod(k, howoftensave) == 0
        xseq((k/howoftensave)+1,:) = x';
    end
    t = t + tstep_calc;
end


%%%%% dynamics of neural activity %%%%%
function diff_x = SWMdend1_fun1(t,x)

global t_stim_start t_stim_end t_final
global m n
global ff re str_ris rid thr_rid
global dend_threshold dend_slope dend_upperbound
global stimulus_signal stimulus_noise_level
global input_ff input_re input_ris input_rid dend_upperbound J_dend

% determine the stimulus at this time (it is changing because of the presumed dynamic noise)
stimulus = stimulus_signal + stimulus_noise_level * randn(1,m);

% whether the stimulus is presented or not
if (t_stim_start <= t) && (t <= t_stim_end)
    current_stimulus = stimulus;
else
    current_stimulus = zeros(1,m);
end

input_ff = ff .* (ones(n,1) * current_stimulus); % n*m matrix
input_re = re * x; % n*1 matrix
input_ris = str_ris * (1/n) * sum(x);
input_rid = max((rid * x - thr_rid), 0); % n*1 matrix

% dendritic saturation, rather than somatic saturation
J_dend = (input_ff + (input_re * ones(1,m))/m) ./ (1 + (input_rid * ones(1,m))); % divisive inhibition on the dendritic input (before thresholding)
dend_output = min(dend_slope .* max(J_dend - dend_threshold, 0), dend_upperbound);
diff_x = -x + max(sum(dend_output,2) - input_ris, 0);
