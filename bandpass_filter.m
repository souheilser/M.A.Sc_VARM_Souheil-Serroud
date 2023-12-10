function filtered_data = butter_bandpass_filter(data, low_cutoff, high_cutoff, Fs, filter_order)
    nyquist_freq = 0.5 * Fs;
    low_norm_cutoff = low_cutoff / nyquist_freq;
    high_norm_cutoff = high_cutoff / nyquist_freq;
    
    [b, a] = butter(filter_order, [low_norm_cutoff, high_norm_cutoff], 'bandpass');
    filtered_data = filtfilt(b, a, data);
end