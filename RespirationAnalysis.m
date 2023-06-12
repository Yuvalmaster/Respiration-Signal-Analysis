close all; clear all; clc; 

%% Load Data
load("B-L08.mat")
Temp = data(:,1);
Respiration = data(:,2);
Fs= 1/(isi/1000);           % Sampling Frequency [Hz]

%% Segments
Normal_resp = Respiration(1:2037);
Normal_temp = Temp(1:2037);

Hyper_resp = Respiration(2038:5082);
Hyper_temp = Temp(2038:5082);

Hypo_resp = Respiration(5082:10064);
Hypo_temp = Temp(5082:10064);

Speech_resp = Respiration(10065:13719);
Speech_temp = Temp(10065:13719);

%% Eupnea Segment
disp('<strong> NORMAL RESPIRATION </strong>')
% Compare Temp Vs. Resp
Time_norm = 0:isi/1000:length(Normal_resp)*isi/1000;             % Time Vector [Sec]
Time_norm = Time_norm(1:end-1); figure;
plot(Time_norm,Normal_resp); hold on; plot(Time_norm,Normal_temp);
title("Normal Respiration: Pneumograph Vs. Temperature sensor signals"); 
xlabel('[sec]'); grid on; legend('Pneumograph - A/D', 'Temperature sensor - Delta Celsius')

% Find peaks and local minimas
Normal_resp_min_vec = Normal_resp(islocalmin(Normal_resp,'MinSeparation',300));    % Normal_resp local minima
Normal_resp_min_loc = find(islocalmin(Normal_resp,'MinSeparation',300)==1);        % Normal_resp local minima
Normal_resp_max_vec = Normal_resp(islocalmax(Normal_resp,'MinSeparation',300));    % Normal_resp local maxima
Normal_resp_max_loc = find(islocalmax(Normal_resp,'MinSeparation',300)==1);        % Normal_resp local maxima

% Plot
figure; plot(Normal_resp); hold on;
scatter(Normal_resp_max_loc,Normal_resp_max_vec)
scatter(Normal_resp_min_loc,Normal_resp_min_vec)
title("Normal Respiration Peak detection"); 
xlabel('[sec]'); ylabel('[A/D]'); grid on;

% Expiration & Inspiration calculation
Expiration_norm_duration = (Normal_resp_min_loc - Normal_resp_max_loc)/Fs;
mean_Exp_norm_dur = mean(Expiration_norm_duration); STD_Exp_dur = std(Expiration_norm_duration);
disp(append(' Mean Normal Expiration Duration: ', num2str(mean_Exp_norm_dur),' [Sec]'))
disp(append(' STD of Normal Expiration Duration: ', num2str(STD_Exp_dur),' [Sec]'))

Inspiration_norm_duration = abs(Normal_resp_min_loc(1:end-1)-Normal_resp_max_loc(2:end))/Fs;
mean_Insp_norm_dur = mean(Inspiration_norm_duration); STD_Insp_dur = std(Inspiration_norm_duration);
disp(append(' Mean Normal Inspiration Duration: ', num2str(mean_Insp_norm_dur),' [Sec]'))
disp(append(' STD of Normal Inspiration Duration: ', num2str(STD_Insp_dur),' [Sec]'))

% Create RR vector over time
Peak_dist_norm=diff(Normal_resp_max_loc)./Fs; %time between peaks in seconds
RR_norm=NaN(1,length(Normal_resp));
RR_norm(1:Normal_resp_max_loc(1))=Peak_dist_norm(1);
RR_norm(Normal_resp_max_loc(end):length(Normal_resp))=Peak_dist_norm(end);
for i=1:length(Normal_resp_max_loc)-1
    RR_norm(Normal_resp_max_loc(i):Normal_resp_max_loc(i+1))=Peak_dist_norm(i);
end

% Plot RR vector over time
figure; plot(Time_norm,60./RR_norm);
title("Respiration Rate Over time: Normal breathing"); 
xlabel('[sec]'); ylabel('[breath/min]'); grid on;

% Breathing cycle over time vs mean
RR_norm_mean = mean(60./RR_norm);
RR_norm_STD = std(60./Peak_dist_norm);
disp(append(' Mean RR of Normal Respiration: ', num2str(RR_norm_mean),' [Breaths/min]'))
disp(append(' STD RR of Normal Respiration: ', num2str(RR_norm_STD),' [Breaths/min]'))
Breathing_cycles_norm = Inspiration_norm_duration + Expiration_norm_duration(2:end);
mean_breathing_cycle_norm = mean(Breathing_cycles_norm);
std_breathing_cycle_norm = std(Breathing_cycles_norm);
disp(append(' Mean Breathing cycle time of Normal Respiration: ', num2str(mean_breathing_cycle_norm),' [Sec]'))
disp(append(' STD Breathing cycle time of Normal Respiration: ', num2str(std_breathing_cycle_norm),' [Sec]'))
disp(' ') ; disp(' -------------------------------------------------------------------------- '); disp(' ') ;

%% Hyperventilation Segment
disp('<strong> Hyperventilation RESPIRATION </strong>')
% Compare Temp Vs. Resp
Time_hyper = 0:isi/1000:length(Hyper_resp)*isi/1000;             % Time Vector [Sec]
Time_hyper = Time_hyper(1:end-1); figure;
plot(Time_hyper,Hyper_resp); hold on; plot(Time_hyper,Hyper_temp);
title("Hyperventilation: Pneumograph Vs. Temperature sensor signals"); 
xlabel('[sec]'); grid on; legend('Pneumograph - A/D', 'Temperature sensor - Delta Celsius')

filtered_hyper_resp = highpass(Hyper_resp,0.5,Fs);               % Filter EMG noises

figure;
plot(Time_hyper,Hyper_resp); hold on; plot(Time_hyper,filtered_hyper_resp) 
title("Filtered vs. Unfiltered signal: Hyperventilation"); 
xlabel('[sec]'); ylabel('[A/D]'); grid on; legend('Unfiltered signal','Filtered signal')

% Peak detection
Hyper_resp_max_vec_b = filtered_hyper_resp(islocalmax(filtered_hyper_resp,'MinSeparation',0));    % hyper_resp local maxima
Hyper_resp_max_loc_b = find(islocalmax(filtered_hyper_resp,'MinSeparation',0)==1);                % hyper_resp local maxima
Hyper_resp_max_loc(1) = Hyper_resp_max_loc_b(1);
Hyper_resp_max_vec(1) = Hyper_resp_max_vec_b(1);

% Remove impossible detections
for i = 1:length(Hyper_resp_max_loc_b)-1
    if abs(Hyper_resp_max_loc_b(i)-Hyper_resp_max_loc_b(i+1)) > 30
        Hyper_resp_max_loc(i) = Hyper_resp_max_loc_b(i);
        Hyper_resp_max_vec(i) = Hyper_resp_max_vec_b(i);
    end
end
Hyper_resp_max_loc = Hyper_resp_max_loc(find(Hyper_resp_max_loc~=0));
Hyper_resp_max_vec = Hyper_resp_max_vec(find(Hyper_resp_max_vec~=0));
[Hyper_resp_max_vec(end+1),Hyper_resp_max_loc(end+1)] = findpeaks(filtered_hyper_resp(2989:end),'MinPeakDistance',55);
Hyper_resp_max_loc(end) = Hyper_resp_max_loc(end) +2989; 

% Plot Filtered signal and scatter peaks detection
figure; plot(filtered_hyper_resp); hold on;
scatter(Hyper_resp_max_loc,Hyper_resp_max_vec)
title("Hyperventilation Peak detection"); 
xlabel('[sec]'); ylabel('[A/D]'); grid on;

% RR over time vector
Peak_dist_hyper=diff(Hyper_resp_max_loc)./Fs; %time between peaks in seconds
RR_hyper=NaN(1,length(filtered_hyper_resp));
RR_hyper(1:Hyper_resp_max_loc(1))=Peak_dist_hyper(1);
RR_hyper(Hyper_resp_max_loc(end):length(filtered_hyper_resp))=Peak_dist_hyper(end);
for i=1:length(Hyper_resp_max_loc)-1
    RR_hyper(Hyper_resp_max_loc(i):Hyper_resp_max_loc(i+1))=Peak_dist_hyper(i);
end


% Plot RR vector over time
figure; plot(Time_hyper,60./RR_hyper);
title("Respiration Rate Over time: Hyperventilation"); 
xlabel('[sec]'); ylabel('[breath/min]'); grid on;

% Breathing cycle over time vs mean
RR_hyper_mean = mean(60./RR_hyper);
RR_hyper_STD = std(60./RR_hyper);
disp(append(' Mean RR of Hyperventilation Respiration: ', num2str(RR_hyper_mean),' [Breaths/min]'))
disp(append(' STD RR of Hyperventilation Respiration: ', num2str(RR_hyper_STD),' [Breaths/min]'))
Breathing_cycles_hyper = Peak_dist_hyper ;
mean_breathing_cycle_hyper = mean(Breathing_cycles_hyper);
std_breathing_cycle_hyper = std(Breathing_cycles_hyper);
disp(append(' Mean Breathing cycle time of Hyperventilation: ', num2str(mean_breathing_cycle_hyper),' [Sec]'))
disp(append(' STD Breathing cycle time of Hyperventilation Respiration: ', num2str(std_breathing_cycle_hyper),' [Sec]'))
disp(' ') ; disp(' -------------------------------------------------------------------------- '); disp(' ') ;

%% Hypoventilation Segment
disp('<strong> Hypoventilation RESPIRATION </strong>')
% Compare Temp Vs. Resp
Time_hypo = 0:isi/1000:length(Hypo_resp)*isi/1000;             % Time Vector [Sec]
Time_hypo = Time_hypo(1:end-1); figure;
plot(Time_hypo,Hypo_resp); hold on; plot(Time_hypo,Hypo_temp);
title("Hypoventilation: Pneumograph Vs. Temperature sensor signals"); 
xlabel('[sec]'); grid on; legend('Pneumograph - A/D', 'Temperature sensor - Delta Celsius')


filtered_hypo_resp = lowpass(Hypo_resp,0.7,Fs);               % Filter EMG noises

[Hypo_resp_max_vec_b,Hypo_resp_max_loc_b] = findpeaks(filtered_hypo_resp,'MinPeakDistance',400);

Hypo_resp_max_loc(1) = Hypo_resp_max_loc_b(1);
Hypo_resp_max_vec(1) = Hypo_resp_max_vec_b(1);


% Remove impossible detections
for i = 1:length(Hypo_resp_max_loc_b)-1
    if abs(Hypo_resp_max_loc_b(i)-Hypo_resp_max_loc_b(i+1)) > 350
        Hypo_resp_max_loc(i) = Hypo_resp_max_loc_b(i);
        Hypo_resp_max_vec(i) = Hypo_resp_max_vec_b(i);
    end
end
Hypo_resp_max_loc = Hypo_resp_max_loc(find(Hypo_resp_max_loc~=0));
Hypo_resp_max_vec = Hypo_resp_max_vec(find(Hypo_resp_max_vec~=0));

[Hypo_resp_max_vec(end+1),Hypo_resp_max_loc(end+1)] = findpeaks(filtered_hypo_resp(4482:end),'MinPeakDistance',400);
Hypo_resp_max_loc(end) = Hypo_resp_max_loc(end) +4482; 

% Plot Filtered signal and scatter peaks detection
figure; plot(filtered_hypo_resp); hold on;
scatter(Hypo_resp_max_loc,Hypo_resp_max_vec)
title("Hypoventilation Peak detection"); 
xlabel('[sec]'); ylabel('[A/D]'); grid on;

% RR over time vector
Peak_dist_hypo = diff(Hypo_resp_max_loc)./Fs; %time between peaks in seconds
RR_hypo=NaN(1,length(filtered_hypo_resp));
RR_hypo(1:Hypo_resp_max_loc(1))=Peak_dist_hypo(1);
RR_hypo(Hypo_resp_max_loc(end):length(filtered_hypo_resp))=Peak_dist_hypo(end);
for i=1:length(Hypo_resp_max_loc)-1
    RR_hypo(Hypo_resp_max_loc(i):Hypo_resp_max_loc(i+1))=Peak_dist_hypo(i);
end


% Plot RR vector over time
figure; plot(Time_hypo,60./RR_hypo);
title("Respiration Rate Over time: Hypoventilation"); 
xlabel('[sec]'); ylabel('[breath/min]'); grid on;

% Breathing cycle over time vs mean
RR_hypo_mean = mean(60./RR_hypo);
RR_hypo_STD = std(60./RR_hypo);
disp(append(' Mean RR of Hypoventilation Respiration: ', num2str(RR_hypo_mean),' [Breaths/min]'))
disp(append(' STD RR of Hypoventilation Respiration: ', num2str(RR_hypo_STD),' [Breaths/min]'))
Breathing_cycles_hypo = Peak_dist_hypo ;
mean_breathing_cycle_hypo = mean(Breathing_cycles_hypo);
std_breathing_cycle_hypo = std(Breathing_cycles_hypo);
disp(append(' Mean Breathing cycle time of Hypoventilation: ', num2str(mean_breathing_cycle_hypo),' [Sec]'))
disp(append(' STD Breathing cycle time of Hypoventilation Respiration: ', num2str(std_breathing_cycle_hypo),' [Sec]'))
disp(' ') ; disp(' -------------------------------------------------------------------------- '); disp(' ') ;
%% Coughing-reading Segment
disp('<strong> Coughing-reading RESPIRATION </strong>')
% Compare Temp Vs. Resp
Time_Speech = 0:isi/1000:length(Speech_temp)*isi/1000;             % Time Vector [Sec]
Time_Speech = Time_Speech(1:end-1); figure;
plot(Time_Speech,Speech_resp); hold on; plot(Time_Speech,Speech_temp);
title("Coughing-reading: Pneumograph Vs. Temperature sensor signals"); 
xlabel('[sec]'); grid on; legend('Pneumograph - A/D', 'Temperature sensor - Delta Celsius')


filtered_speech_temp = highpass(Speech_temp,0.3,Fs);               % Filter EMG noises

figure; plot(Speech_temp); hold on; plot(filtered_speech_temp);
title("Filtered vs. Unfiltered signal: Coughing-reading"); 
xlabel('[sec]'); ylabel('[A/D]'); grid on; legend('Unfiltered signal','Filtered signal')

[speech_temp_vec_b,speech_temp_loc_b] = findpeaks(filtered_speech_temp,'MinPeakDistance',200);

speech_temp_loc(1) = speech_temp_loc_b(1);
speech_temp_vec(1) = speech_temp_vec_b(1);


% Remove impossible detections
for i = 1:length(speech_temp_loc_b)-1
    if abs(speech_temp_loc_b(i)-speech_temp_loc_b(i+1)) > 150
        speech_temp_loc(i) = speech_temp_loc_b(i);
        speech_temp_vec(i) = speech_temp_vec_b(i);
    end
end
speech_temp_loc = speech_temp_loc(find(speech_temp_loc~=0));
speech_temp_vec = speech_temp_vec(find(speech_temp_vec~=0));


% Plot Filtered signal and scatter peaks detection
figure; plot(filtered_speech_temp); hold on;
scatter(speech_temp_loc,speech_temp_vec)
title("Coughing-reading Peak detection"); 
xlabel('[sec]'); ylabel('[A/D]'); grid on;

% RR over time vector
Peak_dist_speech = diff(speech_temp_loc)./Fs; %time between peaks in seconds
RR_speech=NaN(1,length(filtered_speech_temp));
RR_speech(1:speech_temp_loc(1))=Peak_dist_speech(1);
RR_speech(speech_temp_loc(end):length(filtered_speech_temp))=Peak_dist_speech(end);
for i=1:length(speech_temp_loc)-1
    RR_speech(speech_temp_loc(i):speech_temp_loc(i+1))=Peak_dist_speech(i);
end

% Plot RR vector over time
figure; plot(Time_Speech,60./RR_speech);
title("Respiration Rate Over time: Coughing-reading"); 
xlabel('[sec]'); ylabel('[breath/min]'); grid on;

% Breathing cycle over time vs mean
RR_speech_mean = mean(60./RR_speech);
RR_speech_STD = std(60./RR_speech);
disp(append(' Mean RR of Coughing-reading Respiration: ', num2str(RR_speech_mean),' [Breaths/min]'))
disp(append(' STD RR of Coughing-reading Respiration: ', num2str(RR_speech_STD),' [Breaths/min]'))
Breathing_cycles_speech = Peak_dist_speech ;
mean_breathing_cycle_speech = mean(Breathing_cycles_speech);
std_breathing_cycle_speech = std(Breathing_cycles_speech);
disp(append(' Mean Breathing cycle time of Coughing-reading: ', num2str(mean_breathing_cycle_speech),' [Sec]'))
disp(append(' STD Breathing cycle time of Coughing-reading Respiration: ', num2str(std_breathing_cycle_speech),' [Sec]'))
disp(' ') ; disp(' -------------------------------------------------------------------------- '); disp(' ') ;

%% Coherences
figure; wcoherence(Speech_resp,Hypo_resp(1:length(Speech_resp)),Fs)
title("Coherence between Hypoventilation and coughing-reading: pneumograph"); 
figure; wcoherence(Speech_resp(1:length(Hyper_resp)),Hyper_resp,Fs)
title("Coherence between Hyperventilation and coughing-reading: pneumograph"); 

figure; wcoherence(Speech_temp,Hypo_temp(1:length(Speech_temp)),Fs)
title("Coherence between Hypoventilation and coughing-reading: Temperature sensor"); 
figure; wcoherence(Speech_temp(1:length(Hyper_temp)),Hyper_temp,Fs)
title("Coherence between Hyperventilation and coughing-reading: Temperature sensor"); 

%% RR over all segments
RR = [RR_norm,RR_hyper,RR_hypo,RR_speech];
Time = 0:isi/1000:length(RR)*isi/1000; 
Time = Time(1:end-1);


figure; sgtitle("PFT Over time: Normal -> Hyper -> Hypo -> Coughing-reading"); 

subplot(3,1,1); plot(Time,60./RR); hold on;
xline(0,'--',{'Normal'}); xline(2037/Fs,'--',{'Hyperventilation'}); xline(5082/Fs,'--',{'Hyporventilation'}); xline(10064/Fs,'--',{'Coughing-reading'}); 
xlabel('[sec]'); ylabel('[breath/min]'); grid on; title('Respiration Rate')
 
subplot(3,1,2); plot(Time(1:end-1),Respiration(1:13719)); hold on;
xline(0,'--',{'Normal'}); xline(2037/Fs,'--',{'Hyperventilation'}); xline(5082/Fs,'--',{'Hyporventilation'}); xline(10064/Fs,'--',{'Coughing-reading'}); 
xlabel('[sec]'); ylabel('[A/D]'); grid on; title('Pneumography')

subplot(3,1,3); plot(Time(1:end-1),Temp(1:13719)); hold on;
xline(0,'--',{'Normal'}); xline(2037/Fs,'--',{'Hyperventilation'}); xline(5082/Fs,'--',{'Hyporventilation'}); xline(10064/Fs,'--',{'Coughing-reading'}); 
xlabel('[sec]'); ylabel('[Delta Celsius]'); grid on; title('Temperature sensor')