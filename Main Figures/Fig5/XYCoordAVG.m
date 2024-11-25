%%This code plots in a 2D space the actual coordinates of the head of the mouse in every frame recorded by video

%%Load data
filename = 'A047_20191127_2DLC_resnet50_A043Apr6shuffle1_1030000.csv';
[~,~,data] = xlsread(filename) ;
num_data = csvread(filename,3,0);

Pos_head = [num_data(:,2) num_data(:,3) num_data(:,4)];
Pos_back = [num_data(:,5) num_data(:,6) num_data(:,7)];
Pos_tail = [num_data(:,8) num_data(:,9) num_data(:,10)];

%% Save in variables
x_head = Pos_head(:,1); 
y_head = Pos_head(:,2); 
likelihood_head = Pos_head(:,3);

x_back = Pos_back(:,1); 
y_back = Pos_back(:,2); 
likelihood_back = Pos_back(:,3);

x_tail = Pos_tail(:,1); 
y_tail = Pos_tail(:,2); 
likelihood_tail = Pos_tail(:,3);
%% Correction by likelihood weighting

for i = 1:length(x_head) -10

    x_head(i) = sum(x_head(i:i+10).*likelihood_head(i:i+10))/sum(likelihood_head(i:i+10));
    y_head(i) = sum(y_head(i:i+10).*likelihood_head(i:i+10))/sum(likelihood_head(i:i+10));
    
    x_back(i) = sum(x_back(i:i+10).*likelihood_back(i:i+10))/sum(likelihood_back(i:i+10));
    y_back(i) = sum(y_back(i:i+10).*likelihood_back(i:i+10))/sum(likelihood_back(i:i+10));
    
    x_tail(i) = sum(x_tail(i:i+10).*likelihood_tail(i:i+10))/sum(likelihood_tail(i:i+10));
    y_tail(i) = sum(y_tail(i:i+10).*likelihood_tail(i:i+10))/sum(likelihood_tail(i:i+10));
     
end
plot(x_head);
plot(y_head);

plot(x_back);
plot(y_back);

plot(x_tail);
plot(y_tail);
%% Plot spatial positions
% %Downsampling x, y coordinates
%   x_head = [x_head; x_head_2] ;
%   y_head = [y_head; y_head_2];
%  y_head = y_head_POST ;
% % starting temperatures:
figure
t = length(x_head(:,1)):-1:1; 

T = t;
% Make the first frame: 
figure
h_head = scatter(x_head+50*sind(t(1)/5),y_head+50*sind(t(1)/3),30,T*sind(t(1)/2),'filled'); 
% set color axis limits: 
%caxis([500 50])
cb = colorbar; 
colormap(gca,'jet')
ylabel(cb,'time (5 fps)') 
% write the first frame: 
% gif('temperaturedata.gif') <-uncomment to make a gif
% Loop through each subsequent time step: 
for k = 2:length(t) 
     set(h_head,'xdata',x_head+50*sind(t(k)/5),'ydata',y_head+50*sind(t(k)/3),...
        'cdata',T*sind(t(k)/2))
    % pause(0.1) % <-not necessary for making a gif
     % gif  <-uncomment to make a gif
end
xlim([350 1800])
ylim([250 1200])
%% Plot spatial positions
% % starting temperatures:
t = 1:length(x_back(:,1)); 

T = t;
% Make the first frame: 
figure

h_back = scatter(x_back+50*sind(t(1)/5),y_back+50*sind(t(1)/3),15,T*sind(t(1)/2),'filled'); 
% set color axis limits: 
%caxis([20 25]) 
cb = colorbar; 
ylabel(cb,'time (5 fps)') 
% write the first frame: 
% gif('temperaturedata.gif') <-uncomment to make a gif
% Loop through each subsequent time step: 
for k = 2:length(t) 
     set(h_back,'xdata',x_back+50*sind(t(k)/5),'ydata',y_back+50*sind(t(k)/3),...
        'cdata',T*sind(t(k)/2))
    % pause(0.1) % <-not necessary for making a gif
     % gif  <-uncomment to make a gif
end
xlim([300 1700])
ylim([250 1100])
%% Plot spatial positions
% % starting temperatures:
t = 1:length(x_tail(:,1)); 

T = t;
% Make the first frame: 
figure
h_tail = scatter(x_tail+50*sind(t(1)/5),y_tail+50*sind(t(1)/3),15,T*sind(t(1)/2),'filled'); 
% set color axis limits: 
%caxis([20 25]) 
cb = colorbar; 
ylabel(cb,'temperature') 
% write the first frame: 
% gif('temperaturedata.gif') <-uncomment to make a gif
% Loop through each subsequent time step: 
for k = 2:length(t) 
     set(h_tail,'xdata',x_tail+50*sind(t(k)/5),'ydata',y_tail+50*sind(t(k)/3),...
        'cdata',T*sind(t(k)/2))
    % pause(0.1) % <-not necessary for making a gif
     % gif  <-uncomment to make a gif
end

%% Compute displacement

Dist_head = zeros(length(x_head),1);
% Dist_back = zeros(length(x_back),1);
% Dist_tail = zeros(length(x_tail),1);

for i = 2:length(x_head)
    Dist_head(i) = sqrt((x_head(i)-x_head(i-1))^2 + (y_head(i)-y_head(i-1))^2);
%     Dist_back(i) = sqrt((x_back(i)-x_back(i-1))^2 + (y_back(i)-y_back(i-1))^2);
%     Dist_tail(i) = sqrt((x_tail(i)-x_tail(i-1))^2 + (y_tail(i)-y_tail(i-1))^2);
end

Dist_head_noise_idx = find(Dist_head > 20 | Dist_head <0);
Dist_head(Dist_head_noise_idx) = NaN;

% Dist_back_noise_idx = find(Dist_back > 20 | Dist_back <0);
% Dist_back(Dist_back_noise_idx) = NaN;
% 
% Dist_tail_noise_idx = find(Dist_tail > 20 | Dist_tail <0);
% Dist_tail(Dist_tail_noise_idx) = NaN;


figure
plot(Dist_head,'b');
% hold on
% plot(Dist_back,'c');
% hold on
% plot(Dist_tail,'y');

%% Compute velocity
% 
% step_size = 5;

step_size = 20;
%Average velocity with constant acc

velocity_head = zeros(floor(length(Dist_head)/step_size),1);
% velocity_back = zeros(floor(length(Dist_back)/step_size),1);
% velocity_tail = zeros(floor(length(Dist_tail)/step_size),1);

j = 1;

for i = 1:length(velocity_head)
    
    if (j + step_size) < length(Dist_head)
        
    velocity_head(i) = sum(Dist_head(j:j + step_size-1)) / (step_size);
%     velocity_back(i) = sum(Dist_back(j:j + step_size-1)) / (step_size);
%     velocity_tail(i) = sum(Dist_tail(j:j + step_size-1)) / (step_size);
    
    else
    
    velocity_head(i) = sum(Dist_head(j:length(Dist_head)-1)) / (length(Dist_head)-j);
%     velocity_back(i) = sum(Dist_back(j:length(Dist_head)-1)) / (length(Dist_head)-j);
%     velocity_tail(i) = sum(Dist_tail(j:length(Dist_head)-1))/ (length(Dist_head)-j);
    end
    
    j = j + step_size;
    
end

velocity_head_noise_idx = find(velocity_head > 20 | velocity_head <0);
velocity_head(velocity_head_noise_idx) = NaN;

% velocity_back_noise_idx = find(velocity_back > 20 | velocity_back <0);
% velocity_back(velocity_back_noise_idx) = NaN;
% 
% velocity_tail_noise_idx = find(velocity_tail > 20 | velocity_tail <0);
% velocity_tail(velocity_tail_noise_idx) = NaN;

figure
plot(velocity_head,'b');
hold on
plot(velocity_back,'c');
hold on
plot(velocity_tail,'y');
%%
%% Mark time points when the animal was in the nest
Dist_head_in_nest = zeros(length(Dist_head),2);
Dist_head_in_nest(:,1) = Dist_head;

Dist_back_in_nest = zeros(length(Dist_back),2);
Dist_back_in_nest(:,1) = Dist_back;

Dist_tail_in_nest = zeros(length(Dist_tail),2);
Dist_tail_in_nest(:,1) = Dist_tail;

for i = 1: length(Dist_head)
    if x_head(i) > 1165 && x_head(i) < 1500 && y_head(i) > 590 && y_head(i) < 920
        Dist_head_in_nest(i,2) = 1;
    end
end

for i = 1: length(Dist_back)
    if x_back(i) > 1165 && x_back(i) < 1500 && y_back(i) > 590 && y_back(i) < 920
        Dist_back_in_nest(i,2) = 1;
    end
end

for i = 1: length(Dist_tail)
    if x_tail(i) > 1165 && x_tail(i) < 1500 && y_tail(i) > 590 && y_tail(i) < 920
        Dist_tail_in_nest(i,2) = 1;
    end
end

% for i = 1: length(Dist_head)
%     if x_head(i) > 640 && x_head(i) < 820 && y_head(i) > 400 && y_head(i) < 730
%         Dist_head_in_nest(i,2) = 1;
%     end
% end
% 
% for i = 1: length(Dist_back)
%     if x_back(i) > 640 && x_back(i) < 820 && y_back(i) > 400 && y_back(i) < 730
%         Dist_back_in_nest(i,2) = 1;
%     end
% end
% 
% for i = 1: length(Dist_tail)
%     if x_tail(i) > 640 && x_tail(i) < 820 && y_tail(i) > 400 && y_tail(i) < 730
%         Dist_tail_in_nest(i,2) = 1;
%     end
% end


Dist_mouse_in_nest = zeros(length(Dist_head),1);

for i = 1: length(Dist_mouse_in_nest)
    if Dist_head_in_nest(i,2) == 1 && Dist_back_in_nest(i,2) == 1 && Dist_tail_in_nest(i,2) == 1
        Dist_mouse_in_nest(i) = 1;
    end
end

figure
plot(Dist_head_in_nest(:,1),'b')
hold on;
plot(Dist_back_in_nest(:,1),'c')
hold on;
plot(Dist_tail_in_nest(:,1),'y')
hold on;
plot(Dist_mouse_in_nest,'ro')

%% Fit Dist_mouse_in_nest to velocity

Dist_mouse_in_nest_5fps = zeros(floor(length(Dist_mouse_in_nest)/step_size),1);

j = 1;

for i = 1:length(Dist_mouse_in_nest_5fps)
    
    if (j + step_size) < length(Dist_head)
        
    Dist_mouse_in_nest_5fps(i) = sum(Dist_mouse_in_nest(j:j + step_size-1)) / (step_size);
    
    else
    
    Dist_mouse_in_nest_5fps(i) = sum(Dist_mouse_in_nest(j:length(Dist_mouse_in_nest)-1)) / (length(Dist_mouse_in_nest)-j);
   
    end
    
    j = j + step_size;
    
end
figure
plot(velocity_head,'b');
hold on
plot(velocity_back,'c');
hold on
plot(velocity_tail,'y');
hold on
plot(Dist_mouse_in_nest_5fps,'ro')

%%
%Define sleep onset
sleep_onset_optional = zeros(length(velocity_head),1);
for i = 1: length(velocity_head)
    if i + 717 < length(velocity_head)
    if mean(velocity_head(i:i + 239)) < 0.07
        sleep_onset_optional(i) = i ;
    end
    end
end
sleep_indices = find(sleep_onset_optional>0);
sleep_onset_index = sleep_indices(1);
sleep_onset = sleep_onset_index/5;
%% Define nesting
k = 1;
nesting_optional = zeros(sleep_onset_index,1);
for i = 1:sleep_onset_index + 101
    if Dist_mouse_in_nest_5fps(i) ==1 && sum(Dist_mouse_in_nest_5fps(i:i+100)) == 101
    nesting_optional(k) = i;
    end
    k = k+1;
end
nesting_onset = find(nesting_optional~=0,1,'first');
nest_onset = nesting_onset/5
%% Divide data to sleep/nesting/awake states

sleep = [velocity_head(sleep_onset_index:end) velocity_back(sleep_onset_index:end) velocity_tail(sleep_onset_index:end)];
nesting = [velocity_head(nesting_onset:sleep_onset_index) velocity_back(nesting_onset:sleep_onset_index) velocity_tail(nesting_onset:sleep_onset_index)];
wake = [velocity_head(1:nesting_onset) velocity_back(1:nesting_onset) velocity_tail(1:nesting_onset)];

%% Divide nesting to different movements

%nesting = [velocity_head(nesting_onset:5015) velocity_back(nesting_onset:5015) velocity_tail(nesting_onset:5015)];
nesting_time = length(nesting)/5;

sorting = find(nesting(:,1) > 1 & nesting(:,1) < 6);
sorting_time = length(sorting)/5;
pullin = find(nesting(:,1) > 0.3 & nesting(:,2) > 0.3 & (nesting(:,3) <  0.5 | nesting(:,3) > 5));
pullin_time = length(pullin)/5;
suffer = find(nesting(:,1) < 0.2 & nesting(:,2) < 0.2); 
suffer_time = length(suffer)/5;