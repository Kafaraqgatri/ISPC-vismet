%% Code for Day 12

% Launch EEGLab before we begin (or use your
% shortcut on the shortcuts bar)

%%  Reminder of presentation order for function presentations
%        Week1_Order               Week2_Order     
%    ______________________    _____________________
%
%    {'Lyndon Rakusen'    }    {'Angela White'     }
%    {'Jake Dahill-Fuchel'}    {'Mengfei Wu'       }
%    {'Cristian Preciado' }    {'Cagatay Cora'     }
%    {'Weslee Hguyen'     }    {'Robert Poston'    }
%    {'Sebastian Ehmann'  }    {'Ivy Proudfoot'    }
%    {'Jessica Schachtner'}    {'Russell Dougherty'}

%% Geeky bonus for today
% make a plot that fits the day
Four20
% Do the same with a bit of randomness and save it
Four20_rng

    
%% FINALLY, today's long overdue topic: Permutation testing

        
  %% Demonstration with correlation
        
        %fabricate data for correlation for demo
        num_observations = 100;  % mall-to-modest sized sample
        x = randn(1,num_observations)';  % random normal distribution of values (transpose to make column)
        y = x + 9*(rand(1,length(x))');  % by using 9, this gives a population correlation of about .35, the crud factor size  (transpose to make column)
        data = [x y];
        figure(1); scatter(data(:,1),data(:,2));
        true_corr=corr(data(:,1),data(:,2));
        title(['Correlation with fabricated data, ' num2str(num_observations) ' observations, r=' num2str(true_corr)]);
        xlabel('x values');
        ylabel('y values');

        %% Permutations to create distribution of r values under null hypothesis 
        
        % specify permutation parameters
        num_iters = 100;

        % create output vector 1 by number of iterations
        r = nan(1,num_iters);

        % loop through permutations
        for n_iter = 1:num_iters
            data_shuffle = data; % makes a copy on each iteration
            perm_order = randperm(size(data,1));  % random permuted order of indices 1 to number of rows
            data_shuffle(:,2)= data(perm_order,2); % shuffle second column with this random perumted order
            r(n_iter) = corr(data_shuffle(:,1),data_shuffle(:,2)); % store correlation for this iteration
        end

        % find cutoffs
        r_sort = sort(r);

        % p<.05 for positive r
        r_pos_critical = r_sort(round(length(r_sort)*.975));

        % p<.05 for negative r
        r_neg_critical = r_sort(round(length(r_sort)*.025));

        % plot the results
        figure(2); hist(r_sort);
        y_max = ylim; y_max = y_max(2); 
        line([r_pos_critical r_pos_critical],[0 y_max],'color','red');
        line([r_neg_critical r_neg_critical],[0 y_max],'color','red');
        line([true_corr true_corr],[0 y_max],'color','green');
        title(['Fabricated data, ' num2str(num_iters) ' permutations, ' num2str(num_observations) ' observations, r=' num2str(true_corr)]);

        %% Converting to Z scores
        
        % p-value
        pval = 0.025;  % one tailed since we have two tails
        
        % convert p-value to Z value
        z_critical = abs(norminv(pval));
        
        % compute mean and standard deviation maps
        mean_h0 = mean(r);
        std_h0  = std(r);

        % now Z-score observed correlation
        z_observed = (true_corr-mean_h0) / std_h0;

        % Display whether it greater than critical value
        % Location to display
        x_pos = xlim; x_pos=x_pos(1); % min x
        y_pos = ylim; y_pos=y_pos(2); % max y

        if abs(z_observed)>z_critical
            text(x_pos+.1, y_pos-20,['z=' num2str(z_observed,2) ' is significant']);
        else
            text(x_pos+.1, y_pos-20,['z=' num2str(z_observed,2) ' is not significant']);
        end
        
        

        

    %% MINI-ASSIGNMENT #1 -- investigate the parameters
        
        % Alter code above to investigate the impact of the number of iterations [100 1000 10000]
        
        % Using 1000 iterations, change sample size to be smaller or
        % larger [15 50 100 200], generate the correlation, and examine the impact on the significance testing
        
    
    %% ASSIGNMENT #2: Now do with two conditions/groups
     clear all; close all; clc 

    % Fabricate data
    condit1_mean = 5;
    condit2_mean = 10;
    sd1=4;
    sd2=4;
    n1=30;
    n2=30;
    
    % create data -- adding 2*sd*randn (which ranges mostly from +/-2)
    data1= condit1_mean + sd1*randn(n1,1);
    data2= condit2_mean + sd2*randn(n2,1);
    
    obs_diff = mean(data2) - mean(data1);

    % create histograms of the data (hint: hist or histogram)
    figure;
    subplot(1,2,1)
    histogram(data1)
    title('Condition 1')
    xlabel('Value')
    ylabel('Count')

    subplot(1,2,2)
    histogram(data2)
    title('Condition 2')
    xlabel('Value')
    ylabel('Count')

    
    % Now set up null hypothesis and permute
        % If no difference, does not matter which condition data come from
        all_data = [data1; data2];
        n_iter = 5000;
        perm_diffs = zeros(n_iter,1);
        % concatenate matrices vertically
        % Create two Matrices (one for data1, one for data2), to store means from each iteration
        % Do once:
            % derive random order for the concatenated matrix using randperm
            % now compare the mean of each permuted condition
        
        rand_order = randperm(length(all_data));
        perm1 = all_data(rand_order(1:n1));
        perm2 = all_data(rand_order(n1+1:end));
        one_perm_diff = mean(perm2) - mean(perm1);
        
        % Now make a loop and do this n_iter times 
            % derive random order for the concatenated matrix using randperm
            % now compare the mean of each permuted condition
        for i = 1:n_iter
            rand_order = randperm(length(all_data));
            perm1 = all_data(rand_order(1:n1));
            perm2 = all_data(rand_order(n1+1:end));
    
            perm_diffs(i) = mean(perm2) - mean(perm1);
        end
        % and determine whether the true mean difference is greater
        % than expected by chance at p=.05 (two tailed)
            % p-value
            pval = 0.025;  % one tailed since we have two tails
            % convert p-value to Z value
            z_critical = abs(norminv(pval));

            % compute mean and standard deviation of permuted data
            mean_diff = mean(perm_diffs);
            std_diff  = std(perm_diffs);

            % now Z-score observed 
            z_observed = (obs_diff - mean_diff) / std_diff;

            % Is it greater than critical value?
            is_it_significant = abs(z_observed)>z_critical;
            empirical_p = mean(abs(perm_diffs) >= abs(obs_diff));
            
            % plot historgram of permuted diffs vs actual diff
            fprintf('Observed mean difference = %.4f\n', obs_diff);
            fprintf('Mean of permuted diffs   = %.4f\n', mean_diff);
            fprintf('SD of permuted diffs     = %.4f\n', std_diff);
            fprintf('Observed z-score         = %.4f\n', z_observed);
            fprintf('Z-test significant?      = %d\n', is_it_significant);
            fprintf('Empirical permutation p  = %.4f\n', empirical_p);
            

            figure;
            histogram(perm_diffs, 30)
            hold on
            xline(obs_diff, 'r', 'LineWidth', 2)
            title('Permutation Distribution of Mean Differences')
            xlabel('Permuted mean difference')
            ylabel('Count')
            legend('Permuted differences', 'Observed difference')
                
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
        freq_epochs=[EEG.event.type]==100;
        rare_epochs=[EEG.event.type]==101;
        
        EEG_freq_trials = EEG.data(:,:,freq_epochs);
        EEG_rare_trials = EEG.data(:,:,rare_epochs);
        
        % channel and time info
        chan2use = 'Pz';
        P3_timewin = [300 500]; % msec
        chanidx = find(strcmpi({EEG.chanlocs.labels},chan2use));
        P3timeidx = dsearchn(EEG.times',P3_timewin');
        
        % create average of rare for site Pz, all timepoints
        rare_erp = mean(EEG_rare_trials(chanidx,:,:),3);
        rare_erp = squeeze(rare_erp);
        % create average of frequent for site Pz, all timepoints
        freq_erp = mean(EEG_freq_trials(chanidx,:,:),3);
        freq_erp = squeeze(freq_erp);
        
        % plot these ERPs overlaid, with x axis labeled in msec
        figure;
        plot(EEG.times, freq_erp, 'b', 'LineWidth', 2); hold on;
        plot(EEG.times, rare_erp, 'r', 'LineWidth', 2);
        xlabel('Time (ms)');
        ylabel('\muV');
        title(['ERP at ' chan2use]);
        legend('Frequent','Rare');
        xline(P3_timewin(1),'k--');
        xline(P3_timewin(2),'k--');
 
        % from real data, find P3 amp for each condition 
            % hint: max in P3 time window
        freq_P3_amp = max(freq_erp(P3timeidx(1):P3timeidx(2)));
        rare_P3_amp = max(rare_erp(P3timeidx(1):P3timeidx(2)));
        obs_diff = rare_P3_amp - freq_P3_amp;
        % Now set up null hypothesis and permute (using only channel Pz)
            % Hint 1: on each iteration, re-average after shuffling trials
            % Hint 2: on each iteration, take mean on the re-averaged data
        freq_pz = squeeze(EEG_freq_trials(chanidx,:,:));
        rare_pz = squeeze(EEG_rare_trials(chanidx,:,:));
        
        all_pz = [freq_pz rare_pz];
        n_freq = size(freq_pz,2);
        n_rare = size(rare_pz,2);
        n_total = n_freq + n_rare;

        perm_diffs = zeros(n_iter,1);
            % now compare the mean of each permuted condition
        for n = 1:n_iter
            rand_order = randperm(n_total);
    
            perm_freq = all_pz(:,rand_order(1:n_freq));
            perm_rare = all_pz(:,rand_order(n_freq+1:end));
    
            perm_freq_erp = mean(perm_freq,2);
            perm_rare_erp = mean(perm_rare,2);
    
            perm_freq_P3 = max(perm_freq_erp(P3timeidx(1):P3timeidx(2)));
            perm_rare_P3 = max(perm_rare_erp(P3timeidx(1):P3timeidx(2)));
    
            perm_diffs(n) = perm_rare_P3 - perm_freq_P3;
        end

            % do this n_iter times 
            
            % and determine whether the true mean difference is greater
            % than expected by chance at p=.05 (two tailed)
        empirical_p = mean(abs(perm_diffs) >= abs(obs_diff));
        fprintf('Observed P3 difference (rare - freq) = %.4f\n', obs_diff);
        fprintf('Empirical two-tailed p-value = %.4f\n', empirical_p);

        figure;
        histogram(perm_diffs,30); hold on;
        xline(obs_diff,'r','LineWidth',2);
        xlabel('Permuted P3 difference');
        ylabel('Count');
        title('Permutation distribution of P3 amplitude differences');
        legend('Permuted differences','Observed difference');
        % BONUS: Do this for all points and find ranges where two conditions differ accounting for the multiple comparisons
        
        % real ERP difference at each time point
        obs_diff_time = rare_erp - freq_erp;
        perm_diff_time = zeros(length(EEG.times), n_iter);
        perm_max = zeros(n_iter,1);
        perm_min = zeros(n_iter,1);

        for n = 1:n_iter
            rand_order = randperm(n_total);
    
            perm_freq = all_pz(:,rand_order(1:n_freq));
            perm_rare = all_pz(:,rand_order(n_freq+1:end));
    
            perm_freq_erp = mean(perm_freq,2);
            perm_rare_erp = mean(perm_rare,2);
    
            this_diff = perm_rare_erp - perm_freq_erp;
            perm_diff_time(:,n) = this_diff;
    
            perm_max(n) = max(this_diff);
            perm_min(n) = min(this_diff);
        end
        % plot 4 subplots
           % 1 is ERPS
           % 2 is using Z scores from permuted dist to find significance
           % 3 is using extreme values (.95 and .025 at each timepoint) no corrections for multiple comparisons
           % 4 is using extreme values (.95 and .025 at each timepoint) correcting for multiple comparisons
                                       % by taking max and min value across all timepoints on each iteration

        z_crit = abs(norminv(0.025));

        % method 1: Z-score at each timepoint from permuted distribution
        perm_mean_t = mean(perm_diff_time, 2);
        perm_std_t  = std(perm_diff_time, 0, 2);
        z_obs_time  = (obs_diff_time(:) - perm_mean_t) ./ perm_std_t;
        sig_z       = abs(z_obs_time) > z_crit;

        % method 2: per-timepoint percentile thresholds (uncorrected)
        upper_uncorr = prctile(perm_diff_time, 97.5, 2);
        lower_uncorr = prctile(perm_diff_time, 2.5,  2);
        sig_uncorr   = obs_diff_time(:) > upper_uncorr | obs_diff_time(:) < lower_uncorr;

        % method 3: max/min across time (corrected for multiple comparisons)
        upper_corr = prctile(perm_max, 97.5);
        lower_corr = prctile(perm_min, 2.5);
        sig_corr   = obs_diff_time(:) > upper_corr | obs_diff_time(:) < lower_corr;

        figure;
        subplot(221);
        plot(EEG.times, freq_erp, 'b', 'LineWidth', 2); hold on;
        plot(EEG.times, rare_erp, 'r', 'LineWidth', 2);
        xlabel('Time (ms)'); ylabel('\muV');
        title('ERPs'); legend('Frequent','Rare');

        subplot(222);
        plot(EEG.times, obs_diff_time, 'k', 'LineWidth', 1.5); hold on;
        yl = ylim;
        sig_plot = nan(size(obs_diff_time));
        sig_plot(sig_z) = obs_diff_time(sig_z);
        plot(EEG.times, sig_plot, 'r', 'LineWidth', 3);
        xlabel('Time (ms)'); ylabel('Rare - Freq (\muV)');
        title('Significance via Z-score'); ylim(yl);

        subplot(223);
        plot(EEG.times, obs_diff_time, 'k', 'LineWidth', 1.5); hold on;
        plot(EEG.times, upper_uncorr, 'b--');
        plot(EEG.times, lower_uncorr, 'b--');
        sig_plot = nan(size(obs_diff_time));
        sig_plot(sig_uncorr) = obs_diff_time(sig_uncorr);
        plot(EEG.times, sig_plot, 'r', 'LineWidth', 3);
        xlabel('Time (ms)'); ylabel('Rare - Freq (\muV)');
        title('Percentile thresh (uncorrected)');

        subplot(224);
        plot(EEG.times, obs_diff_time, 'k', 'LineWidth', 1.5); hold on;
        yline(upper_corr, 'b--'); yline(lower_corr, 'b--');
        sig_plot = nan(size(obs_diff_time));
        sig_plot(sig_corr) = obs_diff_time(sig_corr);
        plot(EEG.times, sig_plot, 'r', 'LineWidth', 3);
        xlabel('Time (ms)'); ylabel('Rare - Freq (\muV)');
        title('Max/min thresh (MC-corrected)');


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

      % combine trials from both conditions
      all_trials = cat(3, EEG_rare_trials, EEG_freq_trials);
      n_rare   = size(EEG_rare_trials,3);
      n_freq   = size(EEG_freq_trials,3);
      n_total  = n_rare + n_freq;

      n_iter = 1000;
      n_frex = length(tf_rare.frex);
      n_time = length(tf_rare.times);

      perm_tf_diff = zeros(n_iter, n_frex, n_time);
      perm_max     = zeros(n_iter,1);
      perm_min     = zeros(n_iter,1);

      for n = 1:n_iter
          rand_order = randperm(n_total);

          shuf_rare = all_trials(:,:,rand_order(1:n_rare));
          shuf_freq = all_trials(:,:,rand_order(n_rare+1:end));

          tf_shuf_rare = tf_decomp(shuf_rare, EEG.times, EEG.srate, 'chanlocs', EEG.chanlocs, ...
              'lowfreq',2,'highfreq',30,'numfreq',25,'baseline','dB', ...
              'baseline_ms',[-200 0],'logspaced',true, ...
              'cycle_range',[3 10],'save_trials',false);
          tf_shuf_freq = tf_decomp(shuf_freq, EEG.times, EEG.srate, 'chanlocs', EEG.chanlocs, ...
              'lowfreq',2,'highfreq',30,'numfreq',25,'baseline','dB', ...
              'baseline_ms',[-200 0],'logspaced',true, ...
              'cycle_range',[3 10],'save_trials',false);

          this_diff = squeeze(tf_shuf_rare.power(1,:,:) - tf_shuf_freq.power(1,:,:));
          perm_tf_diff(n,:,:) = this_diff;

          perm_max(n) = max(this_diff(:));
          perm_min(n) = min(this_diff(:));

          if mod(n,50)==0, fprintf('Iteration %d/%d\n', n, n_iter); end
      end

      % multiple-comparisons-corrected thresholds from max/min distributions
      upper_thresh = prctile(perm_max, 97.5);
      lower_thresh = prctile(perm_min, 2.5);

      % threshold the observed diff map
      obs_diff_map = squeeze(tf_diff.power(1,:,:));
      thresh_map = obs_diff_map;
      thresh_map(obs_diff_map < upper_thresh & obs_diff_map > lower_thresh) = 0;

      subplot(224);
      contourf(tf_diff.times, tf_diff.frex, thresh_map, 40, 'linecolor','none')
      title('Thresholded Diff (MC-corrected)')
      set(gca,'xlim',times2use,'clim',[-4 6.5])
      colorbar
      xlabel('Time (ms)'), ylabel('Frequency (Hz)')

        


    
    

    