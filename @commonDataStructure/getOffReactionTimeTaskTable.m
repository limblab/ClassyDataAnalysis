function getOffReactionTimeTaskTable(cds,times)
    %this is a method function for the common_data_structure (cds) class, and
    %should be located in a folder '@common_data_structure' with the class
    %definition file and other method files
    %
    %computes the trial variables for the CO task and composes the trial
    %table in the cds using the task variables and the generic trial times
    %passed in from the calling function. This is intended to be called by 
    %the getTrialTable method of the cds class, rather than directly by a
    %user
    
    %get our word timing for changes in the state machine:
    % Isolate the individual word timestamps
    bumpWordBase = hex2dec('50');
    bumpMask=cds.words.word >= (bumpWordBase) & cds.words.word <= (bumpWordBase+5);
    bumpTimes = cds.words.ts(bumpMask)';
    bumpCodes = cds.words.word(bumpMask)';

    wordOTOn = hex2dec('40');
    otMask=bitand(hex2dec('f0'),cds.words.word) == wordOTOn;
    otOnTimes = cds.words.ts( otMask);
    otOnCodes = cds.words.word( otMask);
    
    wordGo = hex2dec('31');
    goCueTime = cds.words.ts(cds.words.word == wordGo);
    
    wordStim=hex2dec('60');
    stimMask=bitand(hex2dec('f0'),cds.words.word) == wordStim;
    stimOnTimes=cds.words.ts( stimMask );
    stimOffTimes=cds.words.ts(stimMask);
    stimCodeList=cds.words.word( stimMask );

    % try to find stime times via sync line, since there is a delay between
    % when the behavior code outputs the stim word and when stimulation
    % occurs
    flag_found_stim_times = 0;
    for i_analog = 1:numel(cds.analog)
        if(any(strcmpi(cds.analog{i_analog}.Properties.VariableNames,'ainp16')))
            flag_found_stim_times = 1;
            
            stim_on = cds.analog{i_analog}.t(find(diff(cds.analog{i_analog}.ainp16 - mean(cds.analog{i_analog}.ainp16) > 3) > 0.5));
            % remove stim times that are close to each other, only store
            % beginning of train
            
            stim_on_mask = find(diff([0;stim_on]) > 0.3); % shift over by one so we get the first stim of each train
            stimOnTimes = stim_on(stim_on_mask);
            stimOffTimes = [stim_on(stim_on_mask(2:end)-1),stim_on(end)];
        end
    end
    
    if(~flag_found_stim_times)
        stimOnTimes=cds.words.ts( stimMask );
    elseif(numel(stimCodeList) ~= numel(stimOnTimes)) % truncate stimCodeList to match stimTimes
        stimCodeTimes=cds.words.ts( stimMask );
        actualStimOnTimes = [];
        actualStimOffTimes = [];
        for i_stim_time = 1:numel(stimCodeTimes)
            % find nearest future actual time, if it's within 200 ms, keep.
            % otherwise discard
            time_diff = stimOnTimes - stimCodeTimes(i_stim_time);
            time_diff(time_diff < 0) = 100000;
            [~,nearest_idx] = min(time_diff);
            
            if(time_diff(nearest_idx) < 0.2)
                actualStimOnTimes(end+1,1) = stimOnTimes(nearest_idx);
                actualStimOffTimes(end+1,1) = stimOffTimes(nearest_idx);
            end
        end
        
        stimOnTimes = actualStimOnTimes;
        stimOffTimes = actualStimOffTimes;
    end
    
    % try to find audio times via sync line, since there is a delay between
    % when the behavior code outputs the audio-stim word and when audio
    % stim occurs
    wordAudio=hex2dec('80');
    audioMask=bitand(hex2dec('f0'),cds.words.word) == wordAudio;
    audioOnTimes=cds.words.ts( audioMask );
    audioOffTimes = cds.words.ts(audioMask);
    audioCodeList=cds.words.word( audioMask );
    
    flag_found_audio_time = 0;
    for i_analog = 1:numel(cds.analog)
        if(any(strcmpi(cds.analog{i_analog}.Properties.VariableNames,'audioout')))
            flag_found_audio_time = 1;
            
            audio_on = cds.analog{i_analog}.t(find(diff(cds.analog{i_analog}.audioout - mean(cds.analog{i_analog}.audioout) > 3) > 0.5));
            % remove stim times that are close to each other, only store
            % beginning of train
            
            audio_on_mask = find(diff([0;audio_on]) > 0.3); % shift over by one so we get the first stim that occurs
            
            audioOnTimes = audio_on(audio_on_mask);
            audioOffTimes = [audio_on(audio_on_mask(2:end)-1);audio_on(end)];
        end
    end
    
    if(~flag_found_audio_time)
        audioOnTimes=cds.words.ts( audioMask );
        audioOffTimes=cds.words.ts(audioMask);
    elseif(numel(audioCodeList) ~= numel(audioOnTimes)) % truncate stimCodeList to match stimTimes
        audioCodeTimes=cds.words.ts( audioMask );
        actualAudioOnTime = [];
        actualAudioOffTime = [];
        for i_audio_time = 1:numel(audioCodeTimes)
            % find nearest future actual time, if it's within 200 ms, keep.
            % otherwise discard
            time_diff = audioOnTimes - audioCodeTimes(i_audio_time);
            time_diff(time_diff < 0) = 100000;
            [~,nearest_idx] = min(time_diff);
            
            if(time_diff(nearest_idx) < 0.2)
                actualAudioOnTime(end+1,1) = audioOnTimes(nearest_idx);
                actualAudioOffTime(end+1,1) = audioOffTimes(nearest_idx);
            end
        end
        
        audioOnTimes = actualAudioOnTime;
        audioOffTimes = actualAudioOffTime;
    end
    
    
    %preallocate our trial variables:
    numTrials=numel(times.number);
    ctrHold=nan(numTrials,1);
    goCueTimes = nan(numTrials,1);
    tgtAngle=nan(numTrials,1);
    moveTime=nan(numTrials,1);
    
    tgtSize=nan(numTrials,1);
    tgtDist = nan(numTrials,1);
    
    cueDuration=nan(numTrials,1);
    stimTrial=false(numTrials,1);
    stimCode=nan(numTrials,1);
    stimOnTime=nan(numTrials,1);
    stimOffTime=nan(numTrials,1);
    
    audioTrial=false(numTrials,1);
    audioCode=nan(numTrials,1);
    audioOnTime=nan(numTrials,1);
    audioOffTime=nan(numTrials,1);
    
    catchTrial=false(numTrials,1);
    
    
    tgtOnTime = nan(numTrials,1);
    
     
    %get the databurst version:
    dbVersion=cds.databursts.db(1,2);
    skipList=[];
    switch dbVersion
        case 1
            % loop thorugh our trials and build our list vectors:
            for i_trial = 1:numTrials
                %find and parse the current databurst:
                idxDB = find(cds.databursts.ts > times.startTime(i_trial) & cds.databursts.ts<times.endTime(i_trial), 1, 'first');
                if isempty(idxDB)
                    skipList=[skipList,i_trial];
                    continue
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
                ctrHold(i_trial) = bytes2float(cds.databursts.db(idxDB,10:13));
                goCueTimes(i_trial) = bytes2float(cds.databursts.db(idxDB,14:17)); % overwritten later
                tgtAngle(i_trial)=bytes2float(cds.databursts.db(idxDB,18:21));
                moveTime(i_trial)=bytes2float(cds.databursts.db(idxDB,22:25));

                tgtSize(i_trial)=bytes2float(cds.databursts.db(idxDB,26:29));
                tgtDist(i_trial) = bytes2float(cds.databursts.db(idxDB,30:33));

                cueDuration(i_trial)=bytes2float(cds.databursts.db(idxDB,34:37));
                stimTrial(i_trial)=cds.databursts.db(idxDB,38);

                audioTrial(i_trial)=cds.databursts.db(idxDB,39);

                catchTrial(i_trial)=cds.databursts.db(idxDB,40);

                
                %now get things that rely only on words and word timing:
                
                % OT on
                idxOT = find(otOnTimes > times.startTime(i_trial) & otOnTimes < times.endTime(i_trial), 1, 'first');
                if isempty(idxOT)
                    tgtOnTime(i_trial) = nan;
                else
                    tgtOnTime(i_trial) = otOnTimes(idxOT);
                end
                

                %Stim code
                idx = find(stimOnTimes > times.startTime(i_trial) & stimOnTimes < times.endTime(i_trial),1,'first');
                if isempty(idx)
                    stimCode(i_trial) = nan;
                else
                    stimCode(i_trial) = bitand(hex2dec('0f'),stimCodeList(idx));%hex2dec('0f') is a bitwise mask for the trailing bit of the word
                    stimOnTime(i_trial)=  stimOnTimes(idx);
                    stimOffTime(i_trial)= stimOffTimes(idx);
                end
                
                %Audio code
                idx = find(audioOnTimes > times.startTime(i_trial) & audioOnTimes < times.endTime(i_trial),1,'first');
                if isempty(idx)
                    stimCode(i_trial) = nan;
                else
                    audioCode(i_trial) = bitand(hex2dec('0f'),audioCodeList(idx));%hex2dec('0f') is a bitwise mask for the trailing bit of the word
                    audioOnTime(i_trial) = audioOnTimes(idx);
                    audioOffTime(i_trial)= audioOffTimes(idx);
                end
            end
            % go cue
            goCueTimes = max([audioOffTime,stimOffTime],[],2,'omitnan');
            
            %build table:
            trialsTable=table(tgtOnTime,goCueTimes,ctrHold,moveTime,cueDuration,...
                                tgtSize,tgtAngle,tgtDist,...
                                stimTrial,stimCode,stimOnTime,stimOffTime,...
                                audioTrial,audioCode,audioOnTime,audioOffTime,catchTrial,...
                                'VariableNames',{'tgtOnTime','goCueTime','ctrHold','moveTime','cueDuration',...
                                'tgtSize','tgtAngle','tgtDist',...
                                'isStimTrial','stimCode','stimOnset','stimOffset',...
                                'isAudioTrial','audioCode','audioOnset','audioOffset','isCatchTrial'});

            trialsTable.Properties.VariableUnits={'s','s','s','s','s',...
                                                    'cm','deg','cm',...
                                                    'bool','int','s','s',...
                                                    'bool','int','s','s','bool'};
            trialsTable.Properties.VariableDescriptions={'outer target onset time','go cue time','center hold time','time to reach target from end of cue','duration of cue',...
                                                            'size of targets','angle of outer target','distance to outer target from center',...
                                                            'was there stimulation','code in the stim word','when did stim start','when did stim stop',...
                                                            'was there audio stim','code in the audio word','when did audio stim start','when did audio stim end','was this a catch trial'};
            
       
        otherwise
            error('getReactionTimeTaskTable:unrecognizedDBVersion',['the trial table code for RT is not implemented for databursts with version#:',num2str(dbVersion)])
    end
    
    trialsTable=[times,trialsTable];
    trialsTable.Properties.Description='Trial table for the offRT task';
    %sanitize trial table by masking off corrupt databursts with nan's:
    mask= ( trialsTable.ctrHold<0           | trialsTable.ctrHold>10000 | ...
            trialsTable.tgtSize<.000001);
    mask(skipList)=1;
    idx=find(mask);
    for j=5:size(trialsTable,2)
        if ~isempty(find(strcmp({'goCueTime','tgtOnTime','bumpTime','tgtID','bumpID'},trialsTable.Properties.VariableNames{j}),1))
            %skip things that are based on the words, not the databurst
            continue
        end
        if islogical(trialsTable{1,j})
            trialsTable{idx,j}=false;
        else
            trialsTable{idx,j}=nan(size(trialsTable{1,j}));
        end
    end
    
    set(cds,'trials',trialsTable)
    evntData=loggingListenerEventData('getCOTaskTable',[]);
    notify(cds,'ranOperation',evntData)
end