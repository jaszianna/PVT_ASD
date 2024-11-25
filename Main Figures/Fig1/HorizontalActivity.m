%% This code were used in every case for computing displacement, called horizontal activity of the mice
%%This code uses excel files created via DeepLabCut ResNet video analysis


%%Load data
filename = 'LB156_20240711_5DLC_resnet50_Midline_anxietyJan12shuffle1_300000.csv';
[~,~,data] = xlsread(filename) ;
num_data = csvread(filename,3,0);

Pos_head = [num_data(:,8) num_data(:,9) num_data(:,10)];
% Pos_back = [num_data(:,5) num_data(:,6) num_data(:,7)];
% Pos_tail = [num_data(:,8) num_data(:,9) num_data(:,10)];
%% Save in variables
x_head = Pos_head(:,1); 
y_head = Pos_head(:,2); 
likelihood_head = Pos_head(:,3);


% x_back = Pos_back(:,1); 
% y_back = Pos_back(:,2); 
% likelihood_back = Pos_back(:,3);
% 
% x_tail = Pos_tail(:,1); 
% y_tail = Pos_tail(:,2); 
% likelihood_tail = Pos_tail(:,3);
%% Correction by likelihood weighting

for i = 1:length(x_head) -10

    x_head(i) = sum(x_head(i:i+10).*likelihood_head(i:i+10))/sum(likelihood_head(i:i+10));
    y_head(i) = sum(y_head(i:i+10).*likelihood_head(i:i+10))/sum(likelihood_head(i:i+10));
    
%     x_back(i) = sum(x_back(i:i+10).*likelihood_back(i:i+10))/sum(likelihood_back(i:i+10));
%     y_back(i) = sum(y_back(i:i+10).*likelihood_back(i:i+10))/sum(likelihood_back(i:i+10));
%     
%     x_tail(i) = sum(x_tail(i:i+10).*likelihood_tail(i:i+10))/sum(likelihood_tail(i:i+10));
%     y_tail(i) = sum(y_tail(i:i+10).*likelihood_tail(i:i+10))/sum(likelihood_tail(i:i+10));
     
end
plot(x_head);
hold on
plot(y_head);
% hold on
% plot(x_back);
% hold on
% plot(y_back);
% hold on
% plot(x_tail);
% hold on
% plot(y_tail);
%% Compute displacement

Dist_head = zeros(length(x_head),1);
% Dist_back = zeros(length(x_back),1);
% Dist_tail = zeros(length(x_tail),1);

for i = 2:length(x_head) - 10
    
       Dist_head(i) = sqrt((x_head(i)-x_head(i-1))^2 + (y_head(i)-y_head(i-1))^2);
       
   if Dist_head(i)>9 && i >10
    
        Dist_head(i)= mean(Dist_head(i-10:i+10));
    
   end
%     Dist_back(i) = sqrt((x_back(i)-x_back(i-1))^2 + (y_back(i)-y_back(i-1))^2);
%     Dist_tail(i) = sqrt((x_tail(i)-x_tail(i-1))^2 + (y_tail(i)-y_tail(i-1))^2);
end

 %Dist_head_noise_idx = find(Dist_head > 30 | Dist_head <0);
 %Dist_head(Dist_head_noise_idx) = NaN;
% 
% Dist_back_noise_idx = find(Dist_back > 20 | Dist_back <0);
% Dist_back(Dist_back_noise_idx) = NaN;
% 
% Dist_tail_noise_idx = find(Dist_tail > 20 | Dist_tail <0);
% Dist_tail(Dist_tail_noise_idx) = NaN;


figure
plot(Dist_head,'b');
hold on
%plot(Dist_back,'c');
% hold on
% plot(Dist_tail,'y');
%% Compute velocity

step_size = 900;
%step_size = 3600;
%step_size = 1200;

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
% velocity_head_noise_idx = find(velocity_head > 20 | velocity_head <0);
% velocity_head(velocity_head_noise_idx) = NaN;
% 
% velocity_back_noise_idx = find(velocity_back > 20 | velocity_back <0);
% velocity_back(velocity_back_noise_idx) = NaN;
% 
% velocity_tail_noise_idx = find(velocity_tail > 20 | velocity_tail <0);
% velocity_tail(velocity_tail_noise_idx) = NaN;

figure
plot(velocity_head,'b');
hold on
% plot(velocity_back,'c');
% hold on
% plot(velocity_tail,'y');