function [dataResampled]=resampleData(data,desiredFreq)
    %[dataResampled]=resampleData(data,desiredFreq)
    %resampleData takes in a time series column matrix 'data' where the 
    %first column is time, and each subsequent column is a column of data
    %samples. desiredFreq is the desired frequency of the resampled data.
    %
    %resampleData assumes that the samples are spaced evenly in time, and 
    %will not produce reliable results if time sampling is jittery. 
    %reampleData estimates the sampling rate of the original data using 
    %the mode of the inter-sample periods and then uses resample to
    %transform the singal into the resampled one. resampleData also
    %resamples the time to make sure there are no extrapolated points.
    
    
    currentFreq=1/mode(diff(data(:,1)));%SF= sample frequency of data
    %upsample the signal:
    %dataD=interp1(data(:,1),data(:,2:end),t_upsamp);
    %upsample the data so that we can use simple decimation rather than
    %interpolation to downsample
    [p,q]=rat(desiredFreq/currentFreq,.0001);
    
    % detrend first because resample assumes endpoints are 0
    a(1,:) = (data(end,2:end)-data(1,2:end))/(data(end,1)-data(1,1));
    a(2,:) = data(1,2:end);
    dataDetrend = zeros(size(data,1),size(data,2)-1);
    for i = 1:size(data,2)-1
        dataDetrend(:,i) = data(:,i+1)-polyval(a(:,i),data(:,1));
    end
    temp=resample(dataDetrend,p,q);
    
    % get rid of extrapolated points
    resamp_vec = ones(size(data,1),1);
    resamp_vec = downsample(upsample(resamp_vec,p),q);
    ty=downsample(upsample(data(:,1),p),q);
    ty=interp1(find(resamp_vec>0),ty(resamp_vec>0),(1:length(ty))');
    extrap_idx = isnan(ty);
    ty(extrap_idx) = [];
    temp(extrap_idx,:) = [];
    
    dataResampled = zeros(size(temp,1),size(temp,2));
    for i=1:size(dataDetrend,2)
        dataResampled(:,i) = temp(:,i)+polyval(a(:,i),ty(:,1));
    end
    
    dataResampled=[ty dataResampled];
end