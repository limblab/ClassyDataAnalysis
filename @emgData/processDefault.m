function processDefault(emg)
    %this is a method function for the emgData class, and should be saved
    %in the @emgData folder
    %
    %emg.process_default()
    %high pass at 10 Hz, rectify, low pass at 25 Hz
    tmp=emg.data;
    
    % find sampling rate
    samprate = 1/mode(diff(tmp.t));

    [blow,alow] = butter(4,20/samprate);
    [bhigh,ahigh] = butter(4,10/samprate,'high');

    EMGIDX = contains(tmp.Properties.VariableNames,'EMG');
    tmp{:,EMGIDX} = filtfilt(blow,alow,abs(filtfilt(bhigh,ahigh,tmp{:,EMGIDX})));
    
    set(emg,'data',tmp)

    evntData=loggingListenerEventData('processDefault',[]);
    notify(emg,'processedDefault',evntData)
end