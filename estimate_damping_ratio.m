function damping_ratio = estimate_damping_ratio(t, response)
    % Find the peaks of the response
    [pks, locs] = findpeaks(response);
    
    % Calculate the time period between successive peaks
    time_diff = diff(t(locs));
    avg_time_period = mean(time_diff);
    
    % Calculate the logarithmic decrement
    log_decrement = log(pks(1:end-1) ./ pks(2:end));
    avg_log_decrement = mean(log_decrement);
    
    % Estimate the damping ratio
    damping_ratio = avg_log_decrement / (2 * pi * sqrt(1 + (avg_log_decrement/(2*pi))^2));
end