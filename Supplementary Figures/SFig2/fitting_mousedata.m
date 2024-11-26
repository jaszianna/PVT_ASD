%% Custom made script written by Orsolya Szalárdy for analyzing power spectral peaks in NREM sleep
% Required data format: separate file for each participant in .dat format,
% where each column corresponds to an electrode and rows correspond to the
% calculated power spectrum at a given frequency, form 0 Hz to the nyquist frequency
%
clear all;
close all;
%% parameters require user input
fPath = 'h:\work\KOKI_mouse_sleep\ÁBRA_20221111\Control EYFP\PRE\'; % set the path ot the input folder (which contains all individual spectum files in .dat)
outfolder = 'h:\work\KOKI_mouse_sleep\ÁBRA_20221111\Control EYFP\PRE\'; % output folder for figures and .xlsx files
addpath('h:\work\KOKI_mouse_sleep\'); % set the forder where the script is located, or add script manually to the path (then, this line is not required)
electrodeidx = [1]; % set the number of electores
electrode_names = {'Fr'}; % set electrode labels (required for the output only)
% set the text you want to write in first raw of the output .xlsx file in
% the next raw
first_1 = {'slope '}; first_2 = {'intercept '}; first_3 = {'max peak freq '}; first_4 = {'max peak amp '}; first_5 = {'max peak area '}; first_6 = {'GOF '};
frequency=0.5:0.25:48; % set frequency range you want to analyze
freqsteps = 0.25; % frequency step in the spectrum files, this should be defined to calculate the length of the matrices
step_correction = 1; % choose whether you want to interpolate the logarithm steps for fitting line (1 = yes; else = no)
frequency_spindle=(9:0.25:18); % set the spindle frequency range used for peak detection
% set the frequency parameters in between you want to performe the line
% fitting; the linear will be fitted using the range f1:f2 and f3:f4
f1 = 2; % frequency for line fitting
f2 = 6; % frequency for line fitting
f3 = 18; % frequency for line fitting
f4 = 48; % frequency for line fitting

%% parameters do not require user input (calculated from the previous section), and initializein variables
List = dir(fullfile(fPath, '*.dat*'));
frequency_intp=([f1:freqsteps:f2 f3:freqsteps:f4]); % frequency range for fitting power low function of the spindle peak measure
GOF = NaN(length(List),length(electrodeidx)); % goodness of fit
SPINDLE=ones(length(List),length(frequency_spindle),numel(electrodeidx))*88; % initialize data matrix
linepars = NaN(length(List),length(electrodeidx)*2);
maxpeak_freq_all = []; % initialize peak frequency matrix
maxpeak_amp_all = []; % initialize peak amplitude matrix
maxpeak_area_all =[]; % initialize peak area data
fs1 = frequency_spindle(1);
fs2 = frequency_spindle(end);
spectra_nrem_mean = [];
p_mpeak_mean = [];
%% calculate peaks and calculate spectral intercept and slope
for loop = 1:length(List) % loop over each participant
    name=([fPath,List(loop).name]);
    spectra_nrem  = textread(name);
    spectra_nrem_mean = [spectra_nrem_mean,spectra_nrem];
    for elec = 1:numel(electrodeidx) % loop over each electrode
        SPINDLE_DIFF = []; % first derivative
        SPINDLE_DIFF2 = []; % second derivative
        
        if electrodeidx(elec) ==0
            ZC = [NaN];
            DC = [NaN];
            ZC_d2 = [NaN];
        elseif  isnan(spectra_nrem(1,elec))==1
            ZC = [NaN];
            DC = [NaN];
            ZC_d2 = [NaN];
        else
            SPINDLE(loop,:,elec)=(spectra_nrem((fs1/freqsteps+1):(fs2/freqsteps+1),electrodeidx(elec))); % extract power for spindle range
            
            MAXPEAK(loop,:,elec)=(spectra_nrem([(f1/freqsteps+1):(f2/freqsteps+1) ((f3/freqsteps+1):(f4/freqsteps+1))],electrodeidx(elec))); % extract power berween f1:f2 and f3:f4
            
            if step_correction == 1 % interpolate log-scale data
                log_freq = log(frequency_intp);
                step = log_freq(end)-log_freq(end-1);
                wholesp = log(MAXPEAK(loop,:,elec));
                stepchange=[log_freq(1):step:log_freq(end)];
                step_interp = interp1(log_freq,wholesp,stepchange,'pchip');
                ab = dsearchn(stepchange',log(f1)); ac = dsearchn(stepchange',log(f2)); % define the linear-fitting region for the interpolated curve
                ad = dsearchn(stepchange',log(f3)); ae = dsearchn(stepchange',log(f4));% define the linear-fitting region for the interpolated curve
                p_mpeak = polyfit(stepchange([(ab:ac) (ad:ae)]),step_interp([(ab:ac) (ad:ae)]),1);
                linefitrange = stepchange([(ab:ac) (ad:ae)]); % the log_frequency range used for line fitting
                linefitdata = step_interp([(ab:ac) (ad:ae)]); % the log_power used for line fitting
            else
                p_mpeak = polyfit(log(frequency_intp),log(MAXPEAK(loop,:,elec)),1);
            end
            % Evaluate regression with R2
            p_mpeak_mean = [p_mpeak_mean;p_mpeak];
            fitted_line = polyval(p_mpeak,linefitrange);
            R = corrcoef(linefitdata,fitted_line);
            R2 = R(1,2).^2;
            GOF(loop,elec)=R2;
            
            figure(loop)
            plot(log(f1:freqsteps:f4),log((spectra_nrem((f1/freqsteps+1):(f4/freqsteps+1),electrodeidx(elec)))),'b');
            hold on
            plot(log(f1:freqsteps:f2),polyval(p_mpeak,log(f1:freqsteps:f2)),'r');
            plot(log(f3:freqsteps:f4),polyval(p_mpeak,log(f3:freqsteps:f4)),'r');
            plot(log(f2:freqsteps:f3),polyval(p_mpeak,log(f2:freqsteps:f3)),'r:');
            
            xticks([log(2) log(4) log(6) log(8) log(10) log(12)  log(14) log(16) log(20) log(30) log(40) log(50)]);
            xticklabels({'2','4','6','8','10','12','14','16','20','30','40','50'});
            ylim([-1 8]);
            xlabel('Frequency (Hz)'), ylabel('power(log) (\muV^2)');
            title(['NREM power spectum and fitted linear for ',List(loop).name(1:end-4)]);
            saveas(gcf, [outfolder, 'fig_',List(loop).name(1:end-4),'.tif']);
            hold off
            
            linepars(loop,elec)=p_mpeak(1); % this gives the slope
            linepars(loop,(elec+numel(electrodeidx)))=p_mpeak(2); % this gives the itercept
            
            SP = SPINDLE(loop,:,elec); % the row corresponding to the current patient
            for k=1:length(SP)-2
                y = [SP(k) SP(k+1) SP(k+2)];
                x = [frequency_spindle(k) frequency_spindle(k+1) frequency_spindle(k+2)];
                p = polyfit(x,y,2); % quadratic fitting using point n-1, n, n+1
                p1 = [2*p(1) 1*p(2)];  % this gives the derivation of the quadratic function
                derivative = polyval(p1,(x(2))); % calculate the derivative in point n
                if k == 1
                    derivative_01 = polyval(p1,(x(1))); %the first element
                elseif k == length(SP)-2
                    derivative_end = polyval(p1,(x(3))); % the last element
                end
                SPINDLE_DIFF = [SPINDLE_DIFF derivative];
                clear p; clear x; clear y; clear p1; clear p2;
            end
            SPINDLE_DIFF = [derivative_01 SPINDLE_DIFF derivative_end]; % first derivative of each datapoint
            % polynom fitting for the first derivate and calculate the slope of
            % calculate the second derivate)
            for k=1:length(SP)-2
                y = [SPINDLE_DIFF(k) SPINDLE_DIFF(k+1) SPINDLE_DIFF(k+2)];
                x = [frequency_spindle(k) frequency_spindle(k+1) frequency_spindle(k+2)];
                p = polyfit(x,y,2); % quadratic fitting using point n-1, n, n+1
                p1 = [2*p(1) 1*p(2)];  % this gives the derivation of the function
                derivative2 = polyval(p1,(x(2))); % calculate the derivative in point n
                if k == 1
                    derivative2_01 = polyval(p1,(x(1)));% and also in point n-1 for the first
                elseif k == length(SP)-2
                    derivative2_end = polyval(p1,(x(3)));% and also in point n+1 for the last
                end
                
                SPINDLE_DIFF2 = [SPINDLE_DIFF2 derivative2];
                clear p; clear x; clear y; clear p1; clear p2;
            end
            SPINDLE_DIFF2 = [derivative2_01 SPINDLE_DIFF2 derivative2_end]; % second derivative of each datapoint
            
            % find the crossing point of 0
            zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0);
            
            zxidx = zci(SPINDLE_DIFF(1:end-1)); % index of the 0 crossing points of the first derivative
            zxidx_d2 = zci(SPINDLE_DIFF2(1:end-1)); % index of the 0 crossing points of the second derivative
            if isempty(zxidx)
                ZC = [NaN];
                DC = [NaN];
            end
            % find each point where the first derivative is 0 and the second derivative is negative (local maxima)
            a=frequency_spindle;
            for k1 = 1:numel(zxidx)
                idxrng = max([1 zxidx(k1)]):min([zxidx(k1)+1 numel(SPINDLE_DIFF)]);
                xrng = a(idxrng);
                yrng = SPINDLE_DIFF(idxrng);
                yrng2 = SPINDLE_DIFF2(idxrng);
                
                ZC(k1) = interp1( yrng(:), xrng(:), 0, 'linear');
                
                DC(k1) = interp1(xrng(:),yrng2(:),(ZC(k1)),'linear');
            end
            for k2 = 1:numel(zxidx_d2)
                idxrng_2 = max([1 zxidx_d2(k2)]):min([zxidx_d2(k2)+1 numel(SPINDLE_DIFF2)]);
                xrng_2 = a(idxrng_2);
                yrng_2d = SPINDLE_DIFF2(idxrng_2);
                ZC_d2(k2) = interp1( yrng_2d(:), xrng_2(:), 0, 'linear');
            end
        end
        zerocross(elec)={ZC};
        second_deriv(elec)={DC};
        second_deriv_zc(elec)={ZC_d2};
        clear zci;
        clear ZC;
        clear DC;
        clear ZC_d2;
        clear zxidx;
        
        spindle=[];
        maxamplitude=[];
        peak_area = [];
        
        a=cell2mat(second_deriv(1,(elec)));
        idx=find(a<0);
        b=cell2mat(zerocross(1,(elec)));
        if isempty(idx)
            spindle = [NaN];
            peak_area = [NaN];
        end
        % calculate the area under the curve
        for num=1:numel(idx)
            spindle(num)=b(idx(num));
            sd_zc = cell2mat(second_deriv_zc(elec));
            [val,idx_sd]=min(abs(sd_zc-spindle(num)));
            minVal=sd_zc(idx_sd);
            if minVal <= spindle(num)
                if idx_sd == length(sd_zc)
                    val_1 = sd_zc(idx_sd); val_2 = frequency_spindle(end);
                else
                    val_1 = sd_zc(idx_sd); val_2 = sd_zc(idx_sd+1);
                end
            elseif minVal > spindle(num)
                if idx_sd ==1
                    val_1 = frequency_spindle(1); val_2 = sd_zc(idx_sd);
                else
                    val_1 = sd_zc(idx_sd-1); val_2 = sd_zc(idx_sd);
                end
            end
            % calculate area
            x_area = [val_1 val_2];
            x_area_space = linspace(val_1,val_2,20);
            y_area=interp1(frequency_spindle,SPINDLE(loop,:,elec),x_area_space,'pchip');
            y_area_space = log(y_area)+abs(min(log(SPINDLE(loop,:,elec))));
            curve_area = trapz(log(x_area_space), y_area_space);
            
            y_white_1= polyval(p_mpeak,log(val_1)); y_white_2 = polyval(p_mpeak,log(val_2));
            y_white_area = linspace(y_white_1,y_white_2,20);
            y_white_area = y_white_area + abs(min(log(SPINDLE(loop,:,elec))));
            white_area = trapz(log(x_area_space), y_white_area);
            
            peak_area(num) = curve_area-white_area;
            
        end
        spindle_freqs(elec)={spindle};
        peak_area_all(elec)={peak_area};
        clear spindle;
        clear peak_area;
        %         end
        xxx=cell2mat(spindle_freqs(elec)); %frequency
        yyy=cell2mat(peak_area_all(elec)); %area
        
        
        for t = 1:length(xxx) % calculate the maximum amplitude
            maxamplitude(t) = log(interp1(frequency_spindle,SPINDLE(loop,:,elec),(xxx(t)),'pchip'))-(polyval(p_mpeak,log(xxx(t))));
            
        end
        
        peaks = [xxx;maxamplitude;yyy];
        peaks = sortrows(peaks.',2,'descend').';
        
        differences= abs(diff(peaks(1,:)));
        del = find(differences<2);
        if isempty(del)
        else
            for i = 1:length(del)
                if isempty(del)
                else
                    peaks(:,del(1)+1) = [];
                    differences= abs(diff(peaks(1,:)));
                    del = find(differences<0.5);
                end
            end
        end
        
        if size(peaks,2)>=1
            maxpeak_freq = peaks(1,1);
            maxpeak_amp = peaks(2,1);
            maxpeak_area = peaks(3,1);
        else
            maxpeak_freq = [NaN];
            maxpeak_amp = [NaN];
            maxpeak_area = [NaN];
        end
        if size(peaks,2)>=2
            secondmax_freq = peaks(1,2);
            secondmax_amp = peaks(2,2);
            secondmax_area = peaks(3,2);
        else
            secondmax_freq = [NaN];
            secondmax_amp = [NaN];
            secondmax_area = [NaN];
        end
        if size(peaks,2)>=3
            thirdmax_freq = peaks(1,3);
            thirdmax_amp = peaks(2,3);
            thirdmax_area = peaks(3,3);
        else
            thirdmax_freq = [NaN];
            thirdmax_amp = [NaN];
            thirdmax_area = [NaN];
        end
        maxpeak_freq_all(loop,elec) = maxpeak_freq; %
        maxpeak_amp_all(loop,elec) = maxpeak_amp; %
        maxpeak_area_all(loop,elec) = maxpeak_area;
        secondmax_freq_all(loop,elec) = secondmax_freq; %
        secondmax_amp_all(loop,elec) = secondmax_amp; %
        secondmax_area_all(loop,elec) = secondmax_area;
        thirdmax_freq_all(loop,elec) = thirdmax_freq; %
        thirdmax_amp_all(loop,elec) = thirdmax_amp; %
        thirdmax_area_all(loop,elec) = thirdmax_area;
        %         xlswrite(strcat(outfolder,'spindle.xlsx'),exp(spindle_freqs{1,(part)}),cell2mat(patients(i)),num2str(part));
        %         xlswrite(strcat(outfolder,'spindle.xlsx'),(spindle_freqs{1,(elec)}),num2str(elec),num2str(loop));
        %         xlswrite(strcat(outfolder,'area.xlsx'),(peak_area_all{1,(elec)}),num2str(elec),num2str(loop));
        % %         xlswrite(strcat(outfolder,'maxpeak.xlsx'),maxpeak_freq,num2str(elec),num2str(loop));
        clear maxamplitude max_idx maxpeak_freq maxpeak_area;
        
    end
    
    clear spindle_freqs;
    clear peak_area_all;
    clear zerocross;
    clear second_deriv;
    clear second_deriv_zc;
end
% GOF_2 = sqrt(GOF); % calculating mean GOF
% B = GOF_2(~isnan(GOF_2)); % remove NaN values
% M = (tanh(mean(atanh(B))))^2 % ==.5*(log(1+B)-log(1-B)); Fisher Z-transform;; back transformed
% SD = (tanh(std(atanh(B))))^2
% MIN = min(GOF,[],'all') % disp min
% MAX = max(GOF,[],'all') % disp max


fullmatrix = [linepars,maxpeak_freq_all,maxpeak_amp_all,maxpeak_area_all,GOF];
otherpeaks = [secondmax_freq_all,secondmax_amp_all,secondmax_area_all,thirdmax_freq_all,thirdmax_amp_all,thirdmax_area_all];
listname = {};
for i = 1:length(List)
    listname{i} = List(i).name;
end
listname = listname';
first_line = [NaN,strcat(first_1,electrode_names),strcat(first_2,electrode_names),strcat(first_3,electrode_names),strcat(first_4,electrode_names),strcat(first_5,electrode_names),strcat(first_6,electrode_names)];
first_line_otherpeaks = [NaN,repmat(electrode_names,1,6)];
out = [first_line;[listname,num2cell(fullmatrix)]];
out_otherpeaks = [first_line_otherpeaks;[listname,num2cell(otherpeaks)]];
xlswrite(strcat(outfolder,'all_results.xlsx'),out);
xlswrite(strcat(outfolder,'otherpeaks_results.xlsx'),out_otherpeaks);
save(strcat(outfolder,'fullmatrix'),'fullmatrix');

meanspectra1 = nanmean(spectra_nrem_mean,2);
p_mean1 = mean(p_mpeak_mean);
figure(1)
            plot(log(f1:freqsteps:f4),log((meanspectra1((f1/freqsteps+1):(f4/freqsteps+1),electrodeidx(elec)))),'b');
            hold on
            plot(log(f1:freqsteps:f2),polyval(p_mean1,log(f1:freqsteps:f2)),'r');
            plot(log(f3:freqsteps:f4),polyval(p_mean1,log(f3:freqsteps:f4)),'r');
            plot(log(f2:freqsteps:f3),polyval(p_mean1,log(f2:freqsteps:f3)),'r:');
            
            xticks([log(2) log(4) log(6) log(8) log(10) log(12)  log(14) log(16) log(20) log(30) log(40) log(50)]);
            xticklabels({'2','4','6','8','10','12','14','16','20','30','40','50'});
            ylim([-1 8]);
            xlabel('Frequency (Hz)'), ylabel('power(log) (\muV^2)');
            title(['Mean NREM power spectum and fitted linear']);
            saveas(gcf, [outfolder, 'fig_avg.tif']);
 hold off



