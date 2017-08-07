function [fhcal,rotcal,Fy_invert, forceOffsets]=getLabParams(labnum,dateTime,rothandle)
    %wrapper function to contain the logic that selects load cell calibrations
    %and rotation data. Intended for use when converting robot handle load cell
    %data from raw form into room coordinates. Provides load cell calibration
    %matrices, rotation matrices, and a variable used to invert the y axis in
    %the case of data where the load cell was installed upside down. This
    %function is called when preprocessing data as it is loaded into
    %matlab. It is not intended for end user use. Measured force offsets
    %for lab 6 are included as an output, but not handled in calling
    %functions.
    %
    %Each lab with a robot should have a block dedicated to it
    %every time the lab 
    forceOffsets = [];
    if labnum==3 %If lab3 was used for data collection
        % Check date of recording to see if it's before or after the
        % change to force handle mounting.
        if datenum(dateTime) < datenum('5/27/2010')            
            fhcal = [ 0.1019 -3.4543 -0.0527 -3.2162 -0.1124  6.6517; ...
                     -0.1589  5.6843 -0.0913 -5.8614  0.0059  0.1503]';
            rotcal = [0.8540 -0.5202; 0.5202 0.8540];                
            Fy_invert = -1; % old force setup was left hand coordnates.
        elseif datenum(dateTime) < datenum('6/28/2011')
            fhcal = [0.0039 0.0070 -0.0925 -5.7945 -0.1015  5.7592; ...
                    -0.1895 6.6519 -0.0505 -3.3328  0.0687 -3.3321]';
            rotcal = [1 0; 0 1];                
            Fy_invert = 1;
        else
            % Fx,Fy,scaleX,scaleY from ATI calibration file:
            % \\citadel\limblab\Software\ATI FT\Calibration\Lab 3\FT7520.cal
            % fhcal = [Fx;Fy]./[scaleX;scaleY]
            % force_offsets acquired empirically by recording static
            % handle.
            fhcal = [-0.0129 0.0254 -0.1018 -6.2876 -0.1127 6.2163;...
                    -0.2059 7.1801 -0.0804 -3.5910 0.0641 -3.6077]'./1000;
            
            Fy_invert = 1;
            if rothandle
                rotcal = [-1 0; 0 1];  
            else
                rotcal = [1 0; 0 1];  
            end
        end
    elseif labnum==2 %if lab2 was used for data collection
        warning('calc_from_raw_script:Lab2LoadCellCalibration','No one noted what the calibration for the Lab2 robot was, so this processing assumes the same parameters as the original LAB3 values. THE FORCE VALUES RESULTING FROM THIS ANALYSIS MAY BE WRONG!!!!!!!!!!!!!!')
        if datenum(dateTime) < datenum('5/27/2010')            
            fhcal = [ 0.1019 -3.4543 -0.0527 -3.2162 -0.1124  6.6517; ...
                     -0.1589  5.6843 -0.0913 -5.8614  0.0059  0.1503]';
            rotcal = [0.8540 -0.5202; 0.5202 0.8540];                
            Fy_invert = -1; % old force setup was left hand coordnates.
        elseif datenum(out_struct.meta.datetime) < datenum('6/28/2011')
            fhcal = [0.0039 0.0070 -0.0925 -5.7945 -0.1015  5.7592; ...
                    -0.1895 6.6519 -0.0505 -3.3328  0.0687 -3.3321]';
            rotcal = [1 0; 0 1];                
            Fy_invert = 1;
        elseif rothandle
            %included this section for consistency. Old Lab2 files 
            %would never have used a rotated handle
            error('calc_from_raw_script:Lab2RotHandle','the rotate handle option was never used in Lab2. If lab2 has been updated with a loadcell and you are using the handle in a rotated position you need to modify raw2handleforce to reflect this')
        end
    elseif labnum==6 %If lab6 was used for data collection
        if datenum(dateTime) < datenum('5/27/2010')            
        % Fx,Fy,scaleX,scaleY from ATI calibration file:
        % \\citadel\limblab\Software\ATI FT\Calibration\Lab 6\FT16018.cal
        % fhcal = [Fx;Fy]./[scaleX;scaleY]
        % force_offsets acquired empirically by recording static
        % handle.
        if datenum(dateTime) < datenum('07-Mar-2016')
            % Fx,Fy,scaleX,scaleY from ATI calibration file:
            % \\citadel\limblab\Software\ATI FT\Calibration\Lab 6\FT16018.cal
            % fhcal = [Fx;Fy]./[scaleX;scaleY]
            % force_offsets acquired empirically by recording static
            % handle.
            fhcal = [0.02653 0.02045 -0.10720 5.94762 0.20011 -6.12048;...
                    0.15156 -7.60870 0.05471 3.55688 -0.09915 3.44508;...
                    10.01343 0.36172 10.30551 0.39552 10.46860 0.38238;...
                    -0.00146 -0.04159 0.14436 0.02302 -0.14942 0.01492;...
                    -0.16542 -0.00272 0.08192 -0.03109 0.08426 0.03519;...
                    0.00377 -0.09455 0.00105 -0.08402 0.00203 -0.08578]'./1000;
            rotcal = eye(6);
            forceOffsets = [];
        elseif datenum(dateTime) < datenum('17-Jul-2017')
            % Fx,Fy,scaleX,scaleY from ATI calibration file:
            % \\citadel\limblab\Software\ATI FT\Calibration\Lab 6\FT16018.cal
            % fhcal = [Fx;Fy]./[scaleX;scaleY]
            % force_offsets acquired empirically by recording static
            % handle.
            fhcal = [0.02653 0.02045 -0.10720 5.94762 0.20011 -6.12048;...
                    0.15156 -7.60870 0.05471 3.55688 -0.09915 3.44508;...
                    10.01343 0.36172 10.30551 0.39552 10.46860 0.38238;...
                    -0.00146 -0.04159 0.14436 0.02302 -0.14942 0.01492;...
                    -0.16542 -0.00272 0.08192 -0.03109 0.08426 0.03519;...
                    0.00377 -0.09455 0.00105 -0.08402 0.00203 -0.08578]'./1000;
            % rotation of the load cell to match forearm frame
            % (load cell is upside down and slightly rotated)
            theta_off = atan2(3,27); %angle offset of load cell to forearm frame- 3 and 27 are the empirircal measures used to generate the angle
%                         theta_off = 0;
            rotcal = [-cos(theta_off) -sin(theta_off) 0    0             0              0;...
                      -sin(theta_off) cos(theta_off)  0    0             0              0;...
                      0                 0             1    0             0              0;...
                      0                 0             0 -cos(theta_off) -sin(theta_off) 0;...
                      0                 0             0 -sin(theta_off) cos(theta_off)  0;...
                      0                 0             0    0             0              1]'; 
            forceOffsets = [-240.5144  245.3220 -103.0073 -567.6240  332.3762 -591.9336]; %measured 3/17/16
%                         force_offsets = [];
        else
            % replaced lab 6 load cell + handle with lab 3 load cell + handle
            % Fx,Fy,scaleX,scaleY from ATI calibration file:
            % \\citadel\limblab\Software\ATI FT\Calibration\Lab 3\FT7520.cal
            % fhcal = [Fx;Fy]./[scaleX;scaleY]
            % force_offsets acquired empirically by recording static
            % handle.

            % fhcal = [-0.0129 0.0254 -0.1018 -6.2876 -0.1127 6.2163;...
            %         -0.2059 7.1801 -0.0804 -3.5910 0.0641 -3.6077]'./1000;
            error('Load cell from RAW code is untested since load cell was replaced in July 17, 2017')
            fhcal = [-0.06745   0.13235  -0.53124 -32.81043  -0.58791  32.43832;...
                    -1.07432  37.46745  -0.41935 -18.73869   0.33458 -18.82582;...
                    -18.56153   1.24337 -18.54582   0.85789 -18.70268   0.63662;...
                    -0.14634   0.36156 -31.67889   0.77952  32.39412  -0.81438;...
                    36.65668  -1.99599 -19.00259   0.79078 -18.87751   0.31411;...
                    -0.31486  18.88139   0.09343  18.96202  -0.46413  18.94001]'./1000;
            
            Fy_invert = 1;
            if rothandle
                rotcal = diag([-1 1 1 -1 1 1]);  
            else
                rotcal = eye(6);  
            end

        end
        Fy_invert = 1;    
        if rothandle
            error('getLabParams:HandleRotated','Handle rotation not implemented for lab 6')  
        end
    else
        error('getLabParams:BadLab',['lab: ',labnum,' is not configured properly']);
    end
end
