function [open_start, open_end, close_start, close_end] = mark_rs_timepoints(EEG, sname, plot_over)
    if nargin < 2
        sname = '';
        plot_over = false;
    elseif nargin < 3
        plot_over = false;
    end
    
    events = EEG.event;
    Fs = EEG.srate;
    open_start_definite = 0; open_end_definite = 0;
    close_start_definite = 0; close_end_definite = 0;
    open_start_tentative = 0; open_end_tentative = 0;
    close_start_tentative = 0; close_end_tentative = 0;
    
    for ii = 1:length(events)
        if strcmp(events(ii).type, '22')
            if open_start_definite > 0
                close_start_definite = events(ii).latency;
            else
                open_start_definite = events(ii).latency;
            end
        elseif strcmp(events(ii).type, '23')
            if open_end_definite > 0
                close_end_definite = events(ii).latency;
            else
                open_end_definite = events(ii).latency;
            end
        elseif strcmp(events(ii).type, '1001, EyesOpen_Start')
            open_start_tentative = events(ii).latency;
        elseif strcmp(events(ii).type, '1002, EyesOpen_End')
            open_end_tentative = events(ii).latency;
        elseif strcmp(events(ii).type, '1003, EyesClose_Start')
            close_start_tentative = events(ii).latency;
        elseif strcmp(events(ii).type, '1004, EyesClose_End')
            close_end_tentative = events(ii).latency;
        end
    end
        
    if open_start_definite > 0
        open_start = open_start_definite;
    else
        open_start = open_start_tentative;
    end
    
    if open_end_definite > 0 && close_end_definite > 0
        open_end = open_end_definite;
    else
        open_end = open_end_tentative;
    end
    
    if close_start_definite > 0
        close_start = close_start_definite;
    else
        close_start = close_start_tentative;
    end
    
    if close_end_definite > 0
        close_end = close_end_definite;
    else
        close_end = close_end_tentative;
    end

    % handle edge case subjects ** this is important since
    % mark_rs_timepoints doesn't work for early pilot subjects
    if strcmp(sname,'JALL_M_00') || strcmp(sname,'JFIT_F_41') || strcmp(sname,'KROA_F_90') || strcmp(sname,'MDEM_F_42') || strcmp(sname,'TTRA_M_89') || strcmp(sname,'LCOC_F_50')
        open_start = 10*Fs; open_end = 300*Fs;
        close_start = 310*Fs; close_end = 600*Fs;
    end
    
    if plot_over
        axes(gca)
        hold on
        plot([open_start/Fs, open_start/Fs], ylim, '-w', 'LineWidth', 3)
        plot([open_end/Fs, open_end/Fs], ylim, '-w', 'LineWidth', 3)
        plot([close_start/Fs, close_start/Fs], ylim, '-m', 'LineWidth', 3)
        plot([close_end/Fs, close_end/Fs], ylim, '-m', 'LineWidth', 3)
    end
end
