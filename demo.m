%Student¡¯s T Kalman Probability Hypothesis Density Filter
%author:Peng Li
%data:2019/6/10
clear all
clc,close all

% This code is based on the code framework written by Professor BA-NGU VO: 
% http://ba-ngu.vo-au.com/publications.html

% the _common folder needs to be added to the MATLAB path
% parameter Scenario can change the measurement model and tracks
% MC_max_runs determines the Monte Carlo runs

%=============================
% Scenario  
%1 means linear & on outliers---->2 means non-linear & on outliers---->3 means non-linear & outliers
Scenario = 3;  
if Scenario == 1
    if_linear = 1; % linear measurements
elseif Scenario  == 2
    if_linear = 2; % non-linear measurements
elseif Scenario  == 3
    if_linear = 2; % non-linear measurements
end
MC_max_runs = 30;  % Monte Carlo runs

%=================================================
% Parameters of environment
ospa_c = 100;                                     % OSPA parameters
track_time=100;                                   % tracking time
filter_max_number = 7;                            % the number of fitlers
T=1;                                              % detection interval
q_noise=[ 1 ; 1 ];                                 % Process noise
F=[1,T,0,0
    0,1,0,0
    0,0,1,T
    0,0,0,1];                                     % State transition matrix
G=[T^2/2,0
    T,0
    0,T^2/2
    0,T];                                         % Noise transfer matrix
H=[1,0,0,0
   0,0,1,0];                                      % Observation matrix of linear scene
r_noise = [10;10];                                % measurement noise
ourliers_noise = [50;50];                         % measurement noise of outliers
R = diag([r_noise(1,1)^2,r_noise(2,1)^2]);        % parametier R of the filters
lamda=10;                                         % poisson average rate of uniform clutter (per scan)
range_K= [ 0 2000; 0 2000];                       % sensor detection region
pdf_K= lamda/prod(range_K(:,2)-range_K(:,1));     % uniform clutter density
ps=0.99;                                          % Survival probability
pd=0.99;                                          % Detection probability
for i=1:filter_max_number
    plot_data(i).d=zeros(1,track_time);           % OSPA errors
    plot_data(i).n=zeros(1,track_time);           % the number of targets erros
    plot_data(i).t=zeros(1,track_time);           % time cost
end
Xtime=[1 1 21 41
    100 100 80 60];                                             % Survival time of 4 targets
target_position=TARGET_FORM(MC_max_runs,Xtime,track_time,Scenario);      % the ture track of targets
%=======================================================

% start Monte Carlo runs
for MC_runs=1:MC_max_runs
    MC_runs
                                                      
    % get measurements-------> meas1 is Cartesian data for linear filters
    % [1 2 3 4 6], meas2 is BAR data for non-linear filters [5 7]
    
    % linear model
    %====================================
    if if_linear == 1     
        Y_measure=MEASUREMENT(H,r_noise,ourliers_noise,target_position,track_time,pd,lamda,MC_max_runs,Xtime);  
        meas1.K = track_time;
        for i = 1:track_time
            meas1.Z{i} = Y_measure(i).m;
        end
        meas2 = turn_plor_ekf(meas1);
    %=====================================
    
    % non-linear model
    %====================================
    elseif if_linear == 2
        model = gen_model_2(q_noise, R,G,7,lamda,pd);
        if Scenario == 2
               model.outliers_flg = 0; % outliers flag off
        elseif Scenario == 3
               model.outliers_flg = 1; % outliers flag on
        end
        truth = gen_truth_ekf(model); % the the ture track of targets 
        [meas1 meas2]=  gen_meas_ekf(model,truth);
        
        % plot figure when MC_max_runs =1
        if MC_max_runs ==1
            handles= plot_results_ekf(model,truth,meas1,[]);
        end
    %========================================
    end
    
    % plot figure when MC_max_runs =1
    if MC_max_runs==1
        figure(2)
        hold on
        for i=1:track_time
            if size(meas1.Z{i},2)>0
                plot(meas1.Z{i}(1,:),meas1.Z{i}(2,:),'b.'),hold on
                if if_linear == 2
                    axis([-2000 2000 0 2000]);
                elseif if_linear == 1
                    axis([0 2000 0 2000]);
                end
            end
        end
        title('measurement')
        xlabel('x position'),ylabel('y position')
    end
    
    % save target true positions for OSPA
    for i=1:track_time
        x(i).m=[];
        if if_linear == 1
            for j=1:size(Xtime,2)
                if target_position.n(j,i)==1
                    x(i).m=[x(i).m,target_position.m(j).s(:,i-Xtime(1,j)+1)];
                end
            end
        elseif if_linear == 2
            x(i).m=[x(i).m truth.X{i}([1,3],:)];
        end
    end
    
    %==========================================================
    %                    start filtering
    %==========================================================
    
    % Prepare memory
    for huancun_i = 1:filter_max_number
        save_est(huancun_i). d = zeros(1,track_time); % OSPA errors
        save_est(huancun_i). n = zeros(1,track_time); % the number of targets erros
        save_est(huancun_i). t = zeros(1,track_time); % time cost
    end
    
    % filtering the data
    for filter_n= [1 2 3 4 5 6 7] 
        filter_n
        
        % Prepare initial model
        if if_linear == 1
            model = gen_model_1(q_noise, R,G,filter_n,lamda,pd);
        elseif if_linear == 2
            model = gen_model_2(q_noise, R,G,filter_n,lamda,pd);
            if Scenario == 2
               model.outliers_flg = 0; % outliers flag off
            elseif Scenario == 3
               model.outliers_flg = 1; % outliers flag on
            end
        end
        
        % filtering
        if filter_n == 1
            est = run_filter_gms(model,meas1); % GM-PHD
        elseif filter_n == 2
            est = run_filter_tkf(model,meas1); % TKF-PHD
        elseif filter_n == 3
            est = run_filter_j_lmb(model,meas1); % J-LMB
        elseif filter_n == 4
            est = run_filter_cphd(model,meas1); % GM-CPHD
        elseif filter_n == 5
            if if_linear == 2
                est = run_filter_tkf_nonliner(model,meas2); % E-TKF-PHD
            else
                for f5 = 1 : 100
                    est.X{f5} = [];
                end
                est.t = zeros(1,100);
                est.N = zeros(1,100);
            end
        elseif filter_n == 6
            est = run_filter_cphd_tkf(model,meas1); % TKF-CPHD
        elseif filter_n == 7
            est = run_filter_ekf(model,meas2); % EKF-PHD
        end
        save_est(filter_n).t=est.t; % save time cost
        
        % plot figure if MC_max_runs = 1
        if MC_max_runs == 1
            for i = 1:track_time
                figure(1)
                hold on
                fuhao = {'b.','ro','y>','g*','ks','<m','+c'};
                if size(est.X{i},2) > 0
                    plot(est.X{i}(1,:),est.X{i}(3,:),fuhao{filter_n});
                end
            end
        end
        
        % calculate OSPA errors
        for i = 1:track_time
            if size(est.X{i},2)>0
                d(i) = OSPA(est.X{i}([1,3],:),x(i).m,ospa_c,2);
            else
                d(i) = ospa_c;
            end
        end
        
        % save OSPA data
        save_est(filter_n).d=d;
        save_est(filter_n).n=est.N';
        
    end
    
    % save data for plot section
    for p_i = 1:filter_n
        plot_data(p_i).d=plot_data(p_i).d + save_est(p_i).d;
        plot_data(p_i).n=plot_data(p_i).n + save_est(p_i).n;
        plot_data(p_i).t=plot_data(p_i).t + save_est(p_i).t;
    end

end

%==========================================================
%                    plot figures
%==========================================================

for i=1:filter_max_number
    plot_data(i).d = plot_data(i).d / MC_max_runs;
    plot_data(i).n = plot_data(i).n / MC_max_runs;
    plot_data(i).t = plot_data(i).t / MC_max_runs;
end
c = {'-b.','-ro','-y>','-g*','-k+','-m','-c'};
if if_linear == 1
    number = sum(target_position.n);
elseif if_linear == 2
    number = truth.N;
end
%-----------------------------------------------------
figure(3)
plot(number),hold on
for i = [1 4 3 7 2 6 5]
    plot(plot_data(i).n,c{i})
end
legend('True number','GM-PHD','GM-CPHD','J-LMB','EKF-PHD','TKF-PHD','TKF-CPHD','E-TKF-PHD')
xlabel('Time (s)'),ylabel('Number of targets')
axis auto;

figure(4)
hold on
c={'b','r','y','g','k','m','c'};
for i = [1 4 3 7 2 6 5]
    plot(plot_data(i).d(1:track_time),c{i})
end
xlabel('Time (s)'),ylabel('OSPA, p=2, c=100')
legend('GM-PHD','GM-CPHD','J-LMB','EKF-PHD','TKF-PHD','TKF-CPHD','E-TKF-PHD')

figure(5)
hold on
c={'-b.','-ro','-y>','-g*','-k+','-m','-c'};
for i = [1 4 3 7 2 6 5]
    plot(plot_data(i).t(1:track_time),c{i})
end
xlabel('Time (s)'),ylabel('Time cost')
legend('GM-PHD','GM-CPHD','J-LMB','EKF-PHD','TKF-PHD','TKF-CPHD','E-TKF-PHD')

