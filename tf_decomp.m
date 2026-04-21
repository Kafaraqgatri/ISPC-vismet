function tf = tf_decomp(EEG_data,EEG_times,EEG_srate,varargin)
% Purpose: Provide full time-frequency decomposition of epoched EEG data
%
% Usage:
%     tf = tf_decomp(EEG.data,EEG.times,EEG.srate);                            % perform tf decomposition with default values
%     tf = tf_decomp(EEG.data,EEG.times,EEG.srate,'chanlocs',EEG.chanlocs);    % also pass chanlocs to output strcuture
%     tf = tf_decomp(EEG.data,EEG.times,EEG.srate,'chanlocs',EEG.chanlocs, ...
%                 ''baseline','pct', 'baseline_ms', [-200 0]);                 % also specify power to be adjusted using percent of baseline
%
% Required inputs:
%     EEG_data       : channels x times x trials EEG data (e.g., EEG.data from EEGLab)
%     EEG_times      : vector of times in msec (e.g., EEG.times from EEGLab)
%     EEG_srate      : sample rate in Hz (e.g., EEG.srate from EEGLab)
%
% Optional inputs (in any order, a combination of keyword and parameter):
%     'chanlocs'     : channel location and name info in EEG.chanlocs format
%     'lowfreq'      : lower frequency cutoff for tf-analysis (default = 3)
%     'highfreq'     : upper frequency cutoff for tf-analysis (default = 40)
%     'numfreq'      : number of frequency bins for tf-analysis (default = 25)
%     'baseline'     : method of baseline correction for power: 'db' (default)
%                        or 'pct' or 'none'
%     'baseline_ms'  : vector of start and stop time in msec for baseline correction
%     'logspaced'    : set to true for log-spaced frequency bins (default) or
%                        false for linear-spaced frequency bins
%     'cycle_range'  : vector of min and max cycles for wavelets
%                        (default is [3 10]); specify same values for fixed
%                        number of cycles across frequencies (e.g. [3 3])
%     'save_trials'  : set to true to save analytical signal for every trial
%                        or to false to return only mean across trials (default)
%     'silentmode'   : supress output that echos parameters being used (default = false)
%
% Ouput:
%   A structure tf with the following obligatory fields
%     tf.power       : channels X frequency X time matrix of power values (mean across trials)
%     tf.phase       : channels X frequency X time matrix of phase values (mean across trials)
%     tf.times       : vector of times (in msec) of tf data. Downsampled by default
%     tf.frex        : vector of frequencies (in Hz) of tf data.
%     tf.params      : structure with parameters used in decomposition
%
%   Structure tf has the following optional fields
%     If channel location info provided
%       tf.chanlocs  : channel location and naming information in EEG.chanlocs format
%                        (althought this is just passed through to the output
%                        structure, so any structure with channel information is permissible)
%     If trial level data are saved,
%       tf.tf_trials : channels X frequency X time x trials analytical signal (complex numbers, no baseline correction)
%
% Example usage specifying all parameters (using default values!)
%   tf = tf_decomp(EEG.data,EEG.times,EEG.srate,'chanlocs',EEG.chanlocs, ...
%                'lowfreq',3,'highfreq',40,'numfreq',25,'baseline','dB', ...
%                'baseline_ms', [-200 0], 'logspaced',true, ...
%                'cycle_range',[3 10],'save_trials',false);
%
% Handy Suggestion:
%   Output can be used with Mike Cohen's tfviewerx as follows:
%     tfviewerx(tf.times,tf.frex,tf.power,tf.chanlocs);
%
%
% Version 1.0, 08-May-2017, John JB Allen (jallen@email.arizona.edu)


%% Make sure this will run or politely return errors and documentation
% Check that arguments are properly specificied

% First Pre-allocate output variable so that if the procedure aborts, we don't get an error
tf = struct();

% specify list of possible valid keywords for this function
valid_keywords = {'chanlocs','lowfreq','highfreq','numfreq','baseline','baseline_ms','logspaced','cycle_range','save_trials','silentmode'};

% Must specify at least 3 and at most 23 arguments
narginchk(3,23);

% next check that varargin arguments come in pairs and are the proper type
for n_param = 1:2:length(varargin) % first, third, fifth etc should be keywords
    validateattributes(varargin{n_param},{'char'},{'nonempty'},mfilename);
    varargin(n_param) = lower(varargin(n_param));   % make lower case to avoid capitalization issues
    if ~ismember(valid_keywords, varargin(n_param))
        fprintf(2,'%s\n',['Error: ''' varargin{n_param} '''' 'is not a valid keyword.']);   %fprintf to device 2 makes font red (no idea why!)
        fprintf(2,'%s\n','Valid keywords are:');
        fprintf(2,'%s\n','  chanlocs, lowfreq, highfreq, numfreq, baseline, baseline_ms, logspaced, cycle_range, save_trials, silentmode');
        doc tf_decomp;
        return
    end
end

% Next parse parameter string and check variable type
% If it is proper, then assign the variables that the users specified
params_proper = true;
for n_param = 1:2:length(varargin) % every other one is a keyword
    switch varargin{n_param}
        case 'chanlocs'                        % if this keyword is found, then we'll assign the variable and check type
            chanlocs = varargin{n_param + 1};  % parameter following the keyword is assigned
            if ~isstruct(chanlocs)             % check variable type, and if wrong, set boolean and set message for later feedback and abort of function
                params_proper = false;
                helpful_message = ['''chanlocs''' ' must be a structure'];
            end
        case 'lowfreq'
            lowfreq = varargin{n_param + 1};
            if ~isscalar(lowfreq) || ~isnumeric(lowfreq)  % a single value that is a number. If it fails either test, set default
                params_proper = false;
                helpful_message = ['''lowfreq''' ' must be a single number'];
            end
        case 'highfreq'
            highfreq = varargin{n_param + 1};
            if ~isscalar(highfreq) || ~isnumeric(highfreq)
                params_proper = false;
                helpful_message = ['''highfreq''' ' must be a single number'];
            end
        case 'numfreq'
            numfreq = varargin{n_param + 1};
            if ~isscalar(numfreq) || ~isnumeric(numfreq)
                params_proper = false;
                helpful_message = ['''numfreq''' ' must be a single number'];
            end
        case 'baseline'
            baseline = varargin{n_param + 1};
            if ~ischar(baseline)
                params_proper = false;
                helpful_message = ['''baseline''' ' must be a string'];
            end
        case 'baseline_ms'
            baseline_ms = varargin{n_param + 1};
            if ~isvector(baseline_ms) || ~isnumeric(baseline_ms)
                params_proper = false;
                helpful_message = ['''basline_ms''' ' must be a vector of two times in msec'];
            end
        case 'logspaced'
            logspaced = varargin{n_param + 1};
            if ~islogical(logspaced)
                params_proper = false;
                helpful_message = ['''logspaced''' ' must be ''true'' or ''false'''];
            end
        case 'cycle_range'
            cycle_range = varargin{n_param + 1};
            if ~isvector(cycle_range) || ~isnumeric(cycle_range)
                params_proper = false;
                helpful_message = ['''cycle_range''' ' must be a vector of two numbers'];
            end
        case 'save_trials'
            save_trials = varargin{n_param + 1};
            if ~islogical(save_trials)
                params_proper = false;
                helpful_message = ['''save_trials''' ' must be ''true'' or ''false'''];
            end
        case 'silentmode'
            silentmode = varargin{n_param + 1};
            if ~islogical(silentmode)
                params_proper = false;
                helpful_message = ['''silentmode''' ' must be ''true'' or ''false'''];
            end
        otherwise
            % not needed since we validated all possible cases above using validateattributes
    end
    if ~params_proper % we found an improper input
        fprintf(2,'%s\n',['Error: ' helpful_message]);  % do this here so that we only have to do this one time instead of in every 'case' section above
        return
    end
end

%% Assign defaults for non-specified parameters
% Assign defaults if user did not specify
if ~exist('silentmode','var')  % Check for silentmode first (needed for all others below)
    silentmode = false;  % don't display anything since this seems rather redundant!
end
if ~exist('chanlocs','var')              % checking iif it is a variable
    chanlocs = struct;                   % assigning default, and then creating user message
    disp_ns('Assigning default empty value to chanlocs',silentmode);
else
    disp_ns('Using chanlocs for output structure',silentmode);
end
if ~exist('lowfreq','var')
    lowfreq = 3;
    disp_ns('Assigning default low frequency cutoff of 3 Hz',silentmode);
else
    disp_ns(['Using low frequency cutoff of ' num2str(lowfreq) ' Hz'],silentmode);
end
if ~exist('highfreq','var')
    highfreq = 40;
    disp_ns('Assigning default high frequency of 40 Hz',silentmode);
else
    disp_ns(['Using high frequency cutoff of ' num2str(highfreq) ' Hz'],silentmode);
end
if ~exist('numfreq','var')
    numfreq = 25;
    disp_ns('Assigning default number of frequencies to 25',silentmode);
else
    disp_ns(['Using ' num2str(numfreq) ' frequencies'],silentmode);
end
if ~exist('baseline','var')
    baseline = 'dB';
    disp_ns('Assigning default to apply dB correction',silentmode);
else
    disp_ns(['Using baseline correction method: ' baseline],silentmode);
end
if ~exist('baseline_ms','var')
    if EEG_times(1) < 0 && EEG_times(end) >= 0      % if 0 mesec is in the interval and there is pre-event data
        baseline_ms = [max(EEG_times(1),-200) 0] ;  % -200 to zero, unless min EEG time is higher than -200, then min EEG time
    else
        baseline_ms = [EEG_times(1)  EEG_times(end)] ; % But if 0 is not in the interval or no pre-event data, then use whole epoch for baseline
    end
    disp_ns(['Assigning default basline ' mat2str(baseline_ms) ' msec'],silentmode);
else
    disp_ns(['Using ' mat2str(baseline_ms) ' msec for baseline correction'],silentmode);
end
if ~exist('logspaced','var')
    logspaced = true;
    disp_ns('Assigning default to log spacing of frequencies',silentmode);
elseif logspaced == true
    disp_ns('Using log spacing of frequencies',silentmode);
else
    disp_ns('Using linear spacing of frequencies',silentmode);
end
if ~exist('cycle_range','var')
    cycle_range = [3 10];
    disp_ns(['Assigning default cycle range  ' mat2str(cycle_range) ' hz'],silentmode);
else
    disp_ns(['Using a range of ' num2str(cycle_range(1)) ' to ' num2str(cycle_range(2)) ' cycles for wavelet creation'],silentmode);
end
if ~exist('save_trials','var')
    save_trials = false;
    disp_ns('Assigning default value to not save decomposition for individual trials ',silentmode);
elseif save_trials == true
    disp_ns('Using input value to save decomposition for individual trials',silentmode);
else
    disp_ns('Using input value to not save decomposition for individual trials',silentmode);
end

%% Error checking
% EEG data should be 3D Chans by Times by Trials
if sum(logical(size(EEG_data))) ~=3
    fprintf(2,'%s\n','EEG data must be 3 dimensional matrix of channels by times by trials');
    return
end

% times specified for baseline should be within range of EEG data

% low freq should be less than high freq

%% Warnings  (FUTURE FEATURE)

% Not enough data to allow for lowest frequency

% very few freqs given the range of freqs


%% Parameters

% what size are we dealing with!
EEG_nbchans =size(EEG_data,1);
EEG_pnts = size(EEG_data,2);
EEG_trials = size(EEG_data,3);

% wavelet parameters

% frequencies
if logspaced
    frex = logspace(log10(lowfreq),log10(highfreq),numfreq);
else
    frex = linspace(lowfreq,highfreq,numfreq);
end

% number of wavelet cycles
wavecycles = logspace(log10(cycle_range(1)),log10(cycle_range(end)),numfreq);

% time base for wavelents
time = -2:1/EEG_srate:2;
% half that time
half_wave_size = (length(time)-1)/2;

% FFT parameters
nWave = length(time);
nData = EEG_pnts*EEG_trials;
nConv = nWave+nData-1;

% Find time indices
% Convert baseline time into indices; note, even if not doing baseline
% correction, default baseline_ms will be initialized above
baseidx = dsearchn(EEG_times',baseline_ms');


%% Wavelet convolution looping over channels and freqs

% initialize matrix for wavelet family
mwf = nan(numfreq,nWave); % mwf = morlet wavelet family

% initialize matrix to hold trial-level time-freq data
tf.tf_trials = zeros(EEG_nbchans,numfreq,EEG_pnts,EEG_trials);

% initialize matrices for power and phase
tf.power = zeros(EEG_nbchans,numfreq,EEG_pnts);
tf.phase = zeros(EEG_nbchans,numfreq,EEG_pnts);

% Let user know how much work there is to do
disp(['Performing wavelet decomposion for ' num2str(EEG_nbchans) ' channels with ' baseline ' baseline correction']);

for n_chan = 1:EEG_nbchans  % loop on channels
    
    % FFT of data (doesn't change on frequency iteration)
    dataX = fft(reshape(EEG_data(n_chan,:,:),1,[]) ,nConv);  % Reshape data to long form for wavelet convolution
    
    for n_freq = 1:numfreq  % loop on frequencies
        if n_chan == 1  % create family of wavelets only once, use for every channel
            s = wavecycles(n_freq)/(2*pi*frex(n_freq));
            mwf(n_freq,:)  = exp(2*1i*pi*frex(n_freq).*time) .* exp(-time.^2./(2*s^2));
        end
        
        % get wavelet FFT
        waveletX = fft(mwf(n_freq,:),nConv);
        waveletX = waveletX ./ max(waveletX); %normalize (important!)
        
        % run convolution
        as = ifft(waveletX.*dataX,nConv);
        as = as(half_wave_size+1:end-half_wave_size);
        as = reshape(as,EEG_pnts,EEG_trials);
        
        % put this into the tf matrix of complex numbers
        tf.tf_trials(n_chan,n_freq,:,:) = as;
        
        % derive mean power across trials
        tf.power(n_chan,n_freq,:) = mean(abs(as).^2 ,2);
        
        % compute ITPC  (it's the absolute length of the mean Eulerized phase
        % angle across trials: abs, mean, euler, angle)
        tf.phase(n_chan,n_freq,:) = abs(mean(exp(1i*angle(as)),2));
        
    end % frequencies
    
    % baseline correct
    baseline_mean = mean(tf.power(n_chan,:,baseidx(1):baseidx(2)),3);
    
    switch lower(baseline) % make all lower case to avoid capitalization errors
        case 'db'  % decibel
            tf.power(n_chan,:,:) = 10*log10( bsxfun(@rdivide, tf.power(n_chan,:,:), baseline_mean) );
        case 'pct'
            tf.power(n_chan,:,:) = 100 * bsxfun(@rdivide, bsxfun(@minus,tf.power(n_chan,:,:),baseline_mean), baseline_mean);
        case 'none'
            % do nothing!
        otherwise
            if n_chan == 1  % display the warning only once, not every channel!
                disp(['For keyword ''baseline'' the parameter ''' baseline ''' is invalid']);
                disp('Valid options are ''db'' ''pct'' or ''none''');
                disp('No baseline correction applied');
            end
    end  % switch -- for baseline correction
    
    % Progress display
    if mod(n_chan,5) == 0
        fprintf('%s',num2str(n_chan)); % lists channel number every 5 channels
    else
        fprintf('%s','.'); % displays a dot for all other channels
    end
    
end % channels

fprintf('\n'); % after displaying channel progres indicator, provide line feed


%% finalize output structure
tf.times = EEG_times;
tf.frex = frex;
tf.params.lowfreq=lowfreq;
tf.params.highfreq=highfreq;
tf.params.numfreq=numfreq;
tf.params.baseline=baseline;
tf.params.baseline_ms=baseline_ms;
tf.params.logspaced=logspaced;
tf.params.cycle_range=cycle_range;
tf.params.date_run = datestr(now);  % might be useful since if this is later saved as a mat file the date could change

% two optional outputs below
if isfield(chanlocs,'labels') % would be true if we passed in EEGlab chanlocs structure
    tf.chanlocs = chanlocs;    % if not true, tf output won't have chanlocs
end
if ~save_trials  % remove this from output if we are not saving trial-level data
    tf = rmfield(tf,'tf_trials');
end

%% Nested function that displays things unless silentmode is set to true
    function disp_ns(message,silentmode)
        if ~silentmode
            disp(message)
        end
    end % disp_ns

%% TODO
% Finish "warning" cell above
% Add log field to tf structure that indicates any warnings
% Add option to save wavelet family to the tf structure?  (Not sure how useful...)
% Avoid a loop
%  Loop for chans: reshape across all chans and trials  (Attempt 1.0 ran but created garbage, and wasn't that much faster!)
%  Loop for freqs: could take advantage of FFT on 2D matrices to do all freqs at once

%% Plotting sample code
% Plot power
% figure;
% contourf(tf.times,tf.frex,squeeze(tf.power(27,:,:)),40,'linecolor','none'); % one channel (27)
% set(gca,'clim',[-5 5],'ydir','normal','xlim',[-300 1000])  % set y and x axis limits

% Plot phase
% figure; 
% contourf(tf.times,tf.frex,squeeze(tf.phase(27,:,:)),40,'linecolor','none');
% set(gca,'clim',[-90 90],'ydir','normal','xlim',[-300 1000])  % set y and x axis limits

% Use with Mike Cohen's tfviewerx
% tfviewerx(tf.times,tf.frex,tf.power,tf.chanlocs)

end % tf_decomp