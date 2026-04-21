%% Answers for Assignment Day 12

clear all; close all; clc;

%% Assignment #2: Now do with two conditions/groups
        % Fabricate data
        condit1_mean = 5;
        condit2_mean = 10;
        sd1=4;
        sd2=4;
        n1=30;
        n2=30;
        n_iter = 1000;
        
        % create data -- adding 2*sd*randn (which ranges mostly from +/-2)
        data1= condit1_mean*ones(n1,1)+2*sd1*randn(n1,1);
        data2= condit2_mean*ones(n2,1)+2*sd2*randn(n2,1);
        
        % create histograms of the data
        hist([data1 data2]);
        % or uncomment to use newer histogram function
        % histogram(data1,'NumBins',12);
        % hold on
        % histogram(data2,'NumBins',12);
        
        % Now set up null hypothesis and permute
            % If no difference, does not matter which condition data come from
            
            % concatenate matrices vertically
            % Thus first 30 are condition 1, second 30 are condition 2
            permute_data = [data1; data2];
            
            % Create two Matrices (one for data1, one for data2), to store means from each iteration
            permute_mean1 = nan(n_iter,1);
            permute_mean2 = nan(n_iter,1);
            
            for perm_n=1:n_iter  % now make a loop and do this n_iter times
                % derive random order for the concatenated matrix using randperm
                permute_data = permute_data(randperm(size(permute_data,1)));
                
                % now compare the mean of each permuted condition
                permute_mean1(perm_n) = mean(permute_data(1:size(data1,1))); % First 30 entries
                permute_mean2(perm_n) = mean(permute_data(size(data1,1)+1:size(permute_data,1))); % second 30 entries
                
            end
             
            
            % and determine whether the true mean difference is greater
            % than expected by chance at p=.05 (two tailed)
                % p-value
                pval = 0.025;  % one tailed since we have two tails
                % convert p-value to Z value
                z_critical = abs(norminv(pval));

                % compute mean and standard deviation of permuted data
                permute_data_diff = permute_mean1 - permute_mean2;
                mean_permute_diff = mean(permute_data_diff);
                std_permute_diff  = std(permute_data_diff);

                % now Z-score observed 
                z_observed = ((condit1_mean-condit2_mean)-mean_permute_diff)/std_permute_diff;

                % Is it greater than critical value?
                is_it_significant = abs(z_observed)>z_critical;
                
                % plot historgram of permuted diffs vs actual diff
                figure; hist(permute_data_diff);
                permute_data_diff = sort(permute_data_diff);
                diff_pos_critical = permute_data_diff(round(length(permute_data_diff)*.975));
                diff_neg_critical = permute_data_diff(round(length(permute_data_diff)*.025));
                
                
                y_max = ylim; y_max = y_max(2); 
                line([diff_pos_critical diff_pos_critical],[0 y_max],'color','red');
                line([diff_neg_critical diff_neg_critical],[0 y_max],'color','red');
                real_data_diff = condit1_mean - condit2_mean;
                line([real_data_diff real_data_diff],[0 y_max],'color','green');
                title(['Fabricated data, ' num2str(n_iter) ' permutations; Z-observed = ' num2str(z_observed)  '; Z-critical = +/-' num2str(z_critical)]);



%% ASSIGNMENT #3 -- determine whether P300 differs between conditions for one subject
    clear all; close all; clc

        % load sample P300 EEG data

        load 'P300.mat'; % change directory as needed


        % COMMENT:  How to generate and save the seed in case you need to do this later!
        n_iter = 1000;
        s = struct; % define structure
        s=rng;  % make sure it has all the fields
        % Need to do this for each iteration!!!
        % put this in your loop below for permutation iterations
        for n=1:n_iter
            s(n)=rng; 
        end
        % after the loop you will need to save the structure with the random seeds
        save('randseed.mat','s');
        % to rerun with same results, later you would need set the seed on each iteration
        load('randseed.mat'); % loads s
        for n=1:n_iter
            rng(s(n)); % set s as seed for this iteration
        end
        % END COMMENT


        % view the data
        pop_eegplot( EEG, 1, 1, 1); % scrolling channels over time showing event codes
        % next figure was generated with EEGLab GUI and eegh to retrieve
        % the command.  It plots trial by trial data at site Pz, sorted by
        % trial type (freq first, rare second)
        figure; pop_erpimage(EEG,1, [10],[[]],'Pz',10,1,{ '100' '101'},[],'type' ,'yerplabel','\muV','erp','on','cbar','on','topo', { [10] EEG.chanlocs EEG.chaninfo } );

        
        % create data for rare and freq
        % logical arrays for indexing trials
        freq_epochs=[EEG.event.type]==100;
        rare_epochs=[EEG.event.type]==101;
        
        % create subsets of data
        EEG_freq_trials = EEG.data(:,:,freq_epochs);
        EEG_rare_trials = EEG.data(:,:,rare_epochs);
        
        % channel and time info
        chan2use = 'Pz';
        P3_timewin = [300 500]; % msec
        chan_idx = find(strcmpi({EEG.chanlocs.labels},chan2use));
        P3timeidx = dsearchn(EEG.times',P3_timewin');
        
        % create average of rare for site Pz, all timepoints
        rare_ERP = mean(EEG_rare_trials(chan_idx,:,:),3);
        
        % create average of frequent for site Pz, all timepoints
        frequent_ERP = mean(EEG_freq_trials(chan_idx,:,:),3);        
        
        % plot these ERPs overlaid, with x axis labeled in msec
        figure; plot(EEG.times,rare_ERP);
        hold on
        plot(EEG.times,frequent_ERP);
        legend({'Rare' 'Freq'});
        xlabel('msec');
        ylabel('\muV');
        % plot positive down
        ax = gca;
        ax.YDir = 'reverse';
        
        
        % from real data, find P3 amp for each condition 
            % hint: max in P3 time window
        rare_P3_ampl = max(rare_ERP(P3timeidx(1):P3timeidx(2)));
        frequent_P3_ampl = max(frequent_ERP(P3timeidx(1):P3timeidx(2)));
        diff_real_data = rare_P3_ampl - frequent_P3_ampl;     
            
        % Now set up null hypothesis and permute
        perm_P3_diff = nan(n_iter,1);
        replicate_prior_analysis = false; % set to true if you saved the mat file randseed.mat
        
        % generate and save the seed in case you need to do this later!
        if replicate_prior_analysis
            load('randseed.mat'); % loads s -- to be used in loop below
        else
            s = struct; % define structure
            s=rng;  % make sure it has all the fields (will update and expand in the iteration loop)
        end
 
            for perm_n =1:n_iter
                
                % Hint 1: on each iteration, re-average after shuffling trials
                all_Pz_data = squeeze(EEG.data(chan_idx,:,:));
                if replicate_prior_analysis
                    rng(s(perm_n)); % set s as seed for this iteration
                else
                    s(perm_n)=rng; % generate and save seed for this iteration
                end

                all_Pz_data = all_Pz_data(:,randperm(size(all_Pz_data,2)));  % shuffles epoch order
                % Mike's solution -- first x epochs, remaining epochs
                    % perm_rare_ERP = mean(all_Pz_data(:,1:size(EEG_rare_trials,3)),2);  % take first set of epochs (# in rare condition) from shuffled trials and consider "rare"
                    % perm_frequent_ERP = mean(all_Pz_data(:,size(EEG_rare_trials,3)+1:size(all_Pz_data,2)),2); % take all subsequent epochs (# in frequent condition) from shuffled trials and consider "frequent"
                % or this works too (and more elegant) -- keep indexing
                % from original, and pull from shuffled trials
                perm_rare_ERP = mean(all_Pz_data(:,rare_epochs),2);  % take first set of epochs (# in rare condition) from shuffled trials and consider "rare"
                perm_frequent_ERP = mean(all_Pz_data(:,freq_epochs),2); % take all subsequent epochs (# in frequent condition) from shuffled trials and consider "frequent"
                
                % DEBUG plot figure to check...
                % figure; plot(EEG.times,perm_rare_ERP); hold on;
                % plot(EEG.times,perm_frequent_ERP);  legend({'Rare' 'Freq'}); 
                
                % Hint 2: on each iteration, take mean difference on the re-averaged data
                perm_rare_P3_ampl = max(perm_rare_ERP(P3timeidx(1):P3timeidx(2)));
                perm_frequent_P3_ampl = max(perm_frequent_ERP(P3timeidx(1):P3timeidx(2)));
                            
                % now compare the mean of each permuted condition
                perm_P3_diff(perm_n) = perm_rare_P3_ampl - perm_frequent_P3_ampl;
                
            end
        
            % now save the structure with the random seeds if this is not a
            % replication analysis
            if ~replicate_prior_analysis
                save('randseed.mat','s');
            end
            
            % plot diffs in P3 amp
            figure; hist(perm_P3_diff);
            
            % and determine whether the true mean difference is greater
            % than expected by chance at p=.05 (two tailed)
            pval = 0.025;  % one tailed since we have two tails
            % convert p-value to Z value
            z_critical = abs(norminv(pval));

            % compute mean and standard deviation 
            mean_permute_diff = mean(perm_P3_diff);
            std_permute_diff  = std(perm_P3_diff);

            % now Z-score observed 
            z_observed = (diff_real_data-mean_permute_diff)/std_permute_diff;

            % Is it greater than critical value?
            if abs(z_observed)>z_critical
                disp(['The P3 amplitude difference of ' num2str(diff_real_data,2) ' uV is significant at p<.05']);
                display_message = ['Site ' chan2use ' observed Z value of ' num2str(z_observed,2) ' is larger than the threshold of ' num2str(z_critical,3)];
                disp(display_message);
            else
                disp(['The P3 amplitude difference of ' num2str(diff_real_data,2) ' uV is NOT significant at p<.05'])
                display_message = ['Site ' chan2use ' observed Z value of ' num2str(z_observed,2) ' is smaller than the threshold of ' num2str(z_critical,3)]
                disp(display_message);
            end

            % find critical values and plot
            perm_P3_diff = sort(perm_P3_diff);
            diff_pos_critical = perm_P3_diff(round(length(perm_P3_diff)*.975));
            diff_neg_critical = perm_P3_diff(round(length(perm_P3_diff)*.025));
            y_max = ylim; y_max = y_max(2); 
            line([diff_pos_critical diff_pos_critical],[0 y_max],'color','red');
            line([diff_neg_critical diff_neg_critical],[0 y_max],'color','red');
            line([diff_real_data diff_real_data], [0 y_max],'color','green');
            title(display_message);
            

            
%% ASSIGMENT #3 BONUS: Do this for all points and find ranges where two conditions differ accounting for the multiple comparisons

% Do not clear variables from above

 % from real data, find amp for each condition at every timepoint
        % we have rare_ERP and frequent_ERP
        diff_real_data = rare_ERP - frequent_ERP;     
            
        % Now set up null hypothesis and permute
        perm_diff = nan(n_iter,length(diff_real_data));
        replicate_prior_analysis = false; % set to true if you saved the mat file randseed.mat
        
        % generate and save the seed in case you need to do this later!
        if replicate_prior_analysis
            load('randseed.mat'); % loads s -- to be used in loop below
        else
            s = struct; % define structure
            s=rng;  % make sure it has all the fields (will update and expand in the iteration loop)
        end
 
            for perm_n =1:n_iter
                
                % On each iteration, re-average after shuffling trials
                all_Pz_data = squeeze(EEG.data(chan_idx,:,:));
                if replicate_prior_analysis
                    rng(s(perm_n)); % set s as seed for this iteration
                else
                    s(perm_n)=rng; % generate and save seed for this iteration
                end

                all_Pz_data = all_Pz_data(:,randperm(size(all_Pz_data,2)));  % shuffles epoch order
                perm_rare_ERP = mean(all_Pz_data(:,rare_epochs),2);  % take first set of epochs (# in rare condition) from shuffled trials and consider "rare"
                perm_frequent_ERP = mean(all_Pz_data(:,freq_epochs),2); % take all subsequent epochs (# in frequent condition) from shuffled trials and consider "frequent"
                
                % DEBUG plot figure to check...
                % figure; plot(EEG.times,perm_rare_ERP); hold on;
                % plot(EEG.times,perm_frequent_ERP);  legend({'Rare' 'Freq'}); 
                
                % On each iteration, take mean difference on the re-averaged data
                            
                % now compare the mean of each permuted condition
                perm_diff(perm_n,:) = perm_rare_ERP - perm_frequent_ERP;
            end
        
            % now save the structure with the random seeds if this is not a
            % replication analysis
            if ~replicate_prior_analysis
                save('randseed.mat','s');
            end
            
                    
            % and determine whether the true mean difference at any point in time is greater
            % than expected by chance at p=.05 (two tailed)
            pval = 0.025;  % one tailed since we have two tails
            % convert p-value to Z value
            z_critical = abs(norminv(pval));

            % compute mean and standard deviation 
            mean_permute_diff = mean(perm_diff,1);
            std_permute_diff  = std(perm_diff,1);

            % now Z-score observed -- NO CORRECTION FOR MULTIPLE COMPARISONS
            z_observed = (diff_real_data-mean_permute_diff)./std_permute_diff;

            % Figure to show z-score vs extremes of permuted distribution,
            % the latter with and without correction for multiple
            % comparisons

            figure;

            % first the ERPs
            subplot(411); plot(EEG.times,rare_ERP);
            title(['Site ' chan2use ' Rare and Freq ERPs'])
            hold on
            plot(EEG.times,frequent_ERP);
            legend({'Rare' 'Freq'});
            xlabel('msec');
            ylabel('\muV');
            % plot positive down
            ax = gca;
            ax.YDir = 'reverse';

            subplot(412) % z-score method
            plot(EEG.times,z_observed);
            title(['Site ' chan2use ' Z-score method of significance']);
            xlabel('msec');
            ylabel('Z Score');
            line(xlim,[z_critical z_critical],'color','red');
            line(xlim,[-z_critical -z_critical],'color','red');
            % Find those points that exceed the critical values 
            significant_points = find(abs(z_observed) > z_critical);
            % Highlight significant points on the z-score plot
            hold on
            scatter(EEG.times(significant_points), z_observed(significant_points), 'filled', 'MarkerFaceColor', 'red');

            subplot(413) % extremes method
            % find critical values and plot
            for n_sample = 1:size(perm_diff,2)
                perm_diff(:,n_sample) = sort(perm_diff(:,n_sample));
            end
            diff_pos_critical = perm_diff(round(length(perm_P3_diff)*.975),:);
            diff_neg_critical = perm_diff(round(length(perm_P3_diff)*.025),:);
            significant_points = diff_real_data > diff_pos_critical | diff_real_data < diff_neg_critical; 
            % plot difference timeseries
            plot(EEG.times,diff_real_data);
            title(['Site ' chan2use ' Extreme values method of significance']);
            xlabel('msec');
            ylabel('Difference \muV');
             % Highlight significant points on the diff plot
            hold on
            scatter(EEG.times(significant_points), diff_real_data(significant_points), 'filled', 'MarkerFaceColor', 'red');
            
            subplot(414) % extremes method adjusting for multiple comparisons
            % find critical values 
            % find min and max values across timepoints for each iteration
            perm_diff_mc = nan(2,size(perm_diff,1));  % 2 = pos and neg, size across second dimension is time
            perm_diff_mc(1,:) = max(perm_diff,[],2); % largest value across all timepoints on each iteration
            perm_diff_mc(2,:) = min(perm_diff,[],2); % smallest value across all timepoints on each iteration
            diff_pos_critical = perm_diff_mc(1,round(length(perm_diff_mc)*.975),:);
            diff_neg_critical = perm_diff_mc(2,round(length(perm_diff_mc)*.025),:);
            significant_points = diff_real_data > diff_pos_critical | diff_real_data < diff_neg_critical; 
            % plot difference timeseries
            plot(EEG.times,diff_real_data);
            title(['Site ' chan2use ' Extreme values method of significance, adjsted for multiple timepoint comparisons']);
            xlabel('msec');
           ylabel('Difference \muV');
             % Highlight significant points on the diff plot
            hold on
            scatter(EEG.times(significant_points), diff_real_data(significant_points), 'filled', 'MarkerFaceColor', 'red');
            
%% ASSIGNMENT #4 -- Extending to TF data

    % You can peak at Mike's code if you want
    edit permutationTesting.m

    % first do TF decomp of P300 data
      clear all; close all; clc
    % load sample P300 EEG data
      load 'P300.mat'; % change directory as needed
    
    % We'll focus on Pz
    chan2use = 'Pz';
    chan_idx = find(strcmpi({EEG.chanlocs.labels},chan2use));
    
    % Separate Rare and Freq
      % logical arrays for indexing trials
      freq_epochs=[EEG.event.type]==100;
      rare_epochs=[EEG.event.type]==101;
      % EEG trials of these types
      EEG_freq_trials = EEG.data(chan_idx,:,freq_epochs);
      EEG_rare_trials = EEG.data(chan_idx,:,rare_epochs);
        
    % Do TF Decompositions
      tf_rare = tf_decomp(EEG_rare_trials,EEG.times,EEG.srate,'chanlocs',EEG.chanlocs, ...
               'lowfreq',2,'highfreq',30,'numfreq',25,'baseline','dB', ...
               'baseline_ms', [-200 0], 'logspaced',true, ...
               'cycle_range',[3 10],'save_trials',false);
      tf_freq = tf_decomp(EEG_freq_trials,EEG.times,EEG.srate,'chanlocs',EEG.chanlocs, ...
               'lowfreq',2,'highfreq',30,'numfreq',25,'baseline','dB', ...
               'baseline_ms', [-200 0], 'logspaced',true, ...
               'cycle_range',[3 10],'save_trials',false);
      tf_diff = tf_rare; tf_diff.power = tf_rare.power - tf_freq.power;
    
      % Show TF Power for channel Pz
      times2use = [-100 800];
      
      figure; 
      subplot(221);
      contourf(tf_rare.times,tf_rare.frex,squeeze(tf_rare.power(1,:,:)),40,'linecolor','none')
      title('Rare')
      set(gca,'xlim',[times2use],'clim',[-4 6.5])
      colorbar
      xlabel('Time (ms)'), ylabel('Frequency (Hz)')
 
      subplot(222);
      contourf(tf_freq.times,tf_freq.frex,squeeze(tf_freq.power(1,:,:)),40,'linecolor','none')
      title('Freq')
      set(gca,'xlim',[times2use],'clim',[-4 6.5])
      colorbar
      xlabel('Time (ms)'), ylabel('Frequency (Hz)')
  
      subplot(223);
      contourf(tf_diff.times,tf_diff.frex,squeeze(tf_diff.power(1,:,:)),40,'linecolor','none')
      title('Diff')
      set(gca,'xlim',[times2use],'clim',[-4 6.5])
      colorbar
      xlabel('Time (ms)'), ylabel('Frequency (Hz)')
    
    
      % Now you take over .... run 1000 iterations to do tf_decomp on
      % shuffled "rare" and "frequent" trials, obtain the difference, and
      % save that difference in the matrix perm_tf_diff.
      % 
      % Then find for each iteration the max and the min tf power value 
      % of the entire diff tf power matrix
      % 
      % then make a thresholded map of the significant tf power values
      % and place it in subplot 4

      % Setup for iterations
      num_iter = 1000;
      num_freq_trials =sum([EEG.event.type]==100);
      num_rare_trials =sum([EEG.event.type]==101);

      % Make big-beautiful-matrix
      perm_tf_diff = nan(2,size(tf_rare.power)); % 2 is max and min

      for n_iter = 1:num_iters
        permuted_trial_order = randperm(size(EEG.data,3));  
        freq_epochs = permuted_trial_order(1:num_freq_trials);
        rare_epochs = permuted_trial_order(num_freq_trials+1:end);
        EEG_freq_trials = EEG.data(chan_idx,:,freq_epochs);
        EEG_rare_trials = EEG.data(chan_idx,:,rare_epochs);
        
        % Do TF Decompositions
        permuted_tf_rare = tf_decomp(EEG_rare_trials,EEG.times,EEG.srate,'chanlocs',EEG.chanlocs, ...
               'lowfreq',2,'highfreq',30,'numfreq',25,'baseline','dB', ...
               'baseline_ms', [-200 0], 'logspaced',true, ...
               'cycle_range',[3 10],'save_trials',false);
        permuted_tf_freq = tf_decomp(EEG_freq_trials,EEG.times,EEG.srate,'chanlocs',EEG.chanlocs, ...
               'lowfreq',2,'highfreq',30,'numfreq',25,'baseline','dB', ...
               'baseline_ms', [-200 0], 'logspaced',true, ...
               'cycle_range',[3 10],'save_trials',false);
        permuted_tf_diff = permuted_tf_rare; permuted_tf_diff.power = permuted_tf_rare.power - permuted_tf_freq.power;
      end

        
      % STILL WORKING!!!
        
        
        
        


    
    

    