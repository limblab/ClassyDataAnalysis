function getTRTTaskTable(cds,times)
    %this is a method function for the common_data_structure (cds) class, and
    %should be located in a folder '@common_data_structure' with the class
    %definition file and other method files
    %
    %cds.getTRTTaskTable(times)
    % returns no value, instead it populates the trials field
    %of the cds assuming the task is a two-workspace random target task. Takes a single
    %input:times, which is a table with 4 columns: number, startTime,
    %endTime, and result. These times define the start and stop of trials
    %as indicated by the state words for trial start and trial end. the
    %result code will be a character 'R':reward 'A':abort 'F':fail
    %'I':incomplete.
    
    corruptDB=0;
    numTrials = length(times.number);
    wordGo = hex2dec('30');
    goCues =  cds.words.ts(bitand(hex2dec('f0'), cds.words.word) == wordGo);

    wordCTHold = hex2dec('A0');
    ctHoldTimes = cds.words.ts(bitand(hex2dec('f0'), cds.words.word) == wordCTHold);
    
    wordTargHold = hex2dec('a0');
    targMask=bitand(hex2dec('f0'),cds.words.word) == wordTargHold;
    targHoldTimes = cds.words.ts( targMask);

    %check DB version number and run appropriate parsing code
    % DB version 0 has 25 bytes before target position
    db_version=cds.databursts.db(1,2);

    if db_version==0
        hdrSize=25;
        numTgt = (cds.databursts.db(1)-hdrSize)/8;

        targStartList=  nan(numTrials,1);
        goCueList=      nan(numTrials,numTgt);
        numTgts=        numTgt*ones(numTrials,1);
        numAttempted=   nan(numTrials,1);
        xOffsets=       nan(numTrials,1); 
        yOffsets=       nan(numTrials,1);
        tgtSizes=       nan(numTrials,1);
        wsnums=         nan(numTrials,1);
        for trial = 1:numel(times.startTime)
            % Find databurst associated with startTime
            dbidx = find(cds.databursts.ts > times.startTime(trial) & cds.databursts.ts < times.endTime(trial));
            if length(dbidx) > 1
                warning('trt_trial_table: multiple databursts @ t = %.3f, using first:%d',times.startTime(trial),trial);
                dbidx = dbidx(1);
            elseif isempty(dbidx)
                warning('trt_trial_table: no/deleted databurst @ t = %.3f, skipping trial:%d',times.startTime(trial),trial);
                corruptDB=1;
                continue;
            end
            if (cds.databursts.db(dbidx,1)-hdrSize)/8 ~= numTgt
                %catch weird/corrupt databursts with different numbers of targets
                warning('trt_trial_table: Inconsistent number of targets @ t = %.3f, skipping trial:%d',times.startTime(trial),trial);
                corruptDB=1;
                continue;
            end

            % Target on times
            idxTargHold = find(targHoldTimes > times.startTime(trial) & targHoldTimes < times.endTime(trial),1,'first');
            %identify trials with corrupt codes that might end up with extra
            %targets
            if isempty(idxTargHold)
                targStart = NaN;
            else
                targStart = targHoldTimes(idxTargHold);
            end
            
            % Go cues
            idxGo = find(goCues > times.startTime(trial) & goCues < times.endTime(trial));

            %get the codes and times for the go cues
            goCue = nan(1,numTgt);
            if isempty(idxGo)
                tgtsAttempted = 0;
            else
                tgtsAttempted = length(idxGo);
            end
            if tgtsAttempted>1
                goCue(1:tgtsAttempted)=goCues(idxGo);
            end

            %identify trials with corrupt end codes that might end up with extra
            %targets
            if length(idxGo) > numTgt
                warning('trt_trial_table: Inconsistent number of targets @ t = %.3f, skipping trial:%d',times.startTime(trial),trial);
                corruptDB=1;
                continue;
            end
            %find target centers
            ctr=bytes2float(cds.databursts.db(dbidx,hdrSize+1:end));
            % Offsets, target size
            xOffset = bytes2float(cds.databursts.db(dbidx,10:13));
            yOffset = bytes2float(cds.databursts.db(dbidx,14:17));
            tgtSize = bytes2float(cds.databursts.db(dbidx,18:21));
            wsnum = bytes2float(cds.databursts.db(dbidx,22:25));

            % Build arrays
            targStartList(trial,:)=     targStart;          % time of first target onset
            goCueList(trial,:)=         goCue;              % time stamps of go_cue(s)
            numTgts(trial)=             numTgt;             % max number of targets
            numAttempted(trial,:)=      tgtsAttempted;      % ?
            xOffsets(trial)=            xOffset;            % x offset
            yOffsets(trial)=            yOffset;            % y offset
            tgtSizes(trial)=            tgtSize;            % target size
            wsnums(trial)=              wsnum;              % workspace number
            tgtCtrs(trial,:)=           ctr;                %center positions of the targets
        end

        trials=table(targStartList,goCueList,numTgts,numAttempted,xOffsets,yOffsets,tgtSizes,wsnums,tgtCtrs,...
                    'VariableNames',{'targetStartTime','goCueTime','numTgt','numAttempted','xOffset','yOffset','tgtSize','spaceNum','tgtCtr'});
        trials.Properties.VariableUnits={'s','s','int','int','cm','cm','cm','int','cm'};
        trials.Properties.VariableDescriptions={'first target hold time','go cue time','number of targets','number of targets attempted','x offset','y offset','target size','workspace number','target center position'};

    elseif db_version==1
        hdrSize=25;
        numTgt = (cds.databursts.db(1)-hdrSize)/8;

        ctHoldList=     nan(numTrials,1);
        targStartList=  nan(numTrials,1);
        goCueList=      nan(numTrials,numTgt);
        numTgts=        numTgt*ones(numTrials,1);
        numAttempted=   nan(numTrials,1);
        xOffsets=       nan(numTrials,1); 
        yOffsets=       nan(numTrials,1);
        tgtSizes=       nan(numTrials,1);
        wsnums=         nan(numTrials,1);
        for trial = 1:numel(times.startTime)
            % Find databurst associated with startTime
            dbidx = find(cds.databursts.ts > times.startTime(trial) & cds.databursts.ts < times.endTime(trial));
            if length(dbidx) > 1
                warning('trt_trial_table: multiple databursts @ t = %.3f, using first:%d',times.startTime(trial),trial);
                dbidx = dbidx(1);
            elseif isempty(dbidx)
                warning('trt_trial_table: no/deleted databurst @ t = %.3f, skipping trial:%d',times.startTime(trial),trial);
                corruptDB=1;
                continue;
            end
            if (cds.databursts.db(dbidx,1)-hdrSize)/8 ~= numTgt
                %catch weird/corrupt databursts with different numbers of targets
                warning('trt_trial_table: Inconsistent number of targets @ t = %.3f, skipping trial:%d',times.startTime(trial),trial);
                corruptDB=1;
                continue;
            end

            % CT hold start times
            idxCTHold = find(ctHoldTimes > times.startTime(trial) & ctHoldTimes < times.endTime(trial));
            %identify trials with corrupt codes that might end up with extra
            %targets
            if isempty(idxCTHold)
                ctHold = NaN;
            elseif length(idxCTHold) > 1
                warning('trt_trial_table: Multiple center hold @ t = %.3f, skipping trial:%d',times.startTime(trial),trial);
                corruptDB=1;
                continue;
            else
                ctHold = ctHoldTimes(idxInitGo);
            end
            
            % Target on times
            idxInitGo = find(initGoCues > times.startTime(trial) & initGoCues < times.endTime(trial));
            %identify trials with corrupt codes that might end up with extra
            %targets
            if isempty(idxInitGo)
                targStart = NaN;
            elseif length(idxInitGo) > 1
                warning('trt_trial_table: Multiple initial go cues @ t = %.3f, skipping trial:%d',times.startTime(trial),trial);
                corruptDB=1;
                continue;
            else
                targStart = initGoCues(idxInitGo);
            end
            
            % Go cues
            idxGo = find(goCues > times.startTime(trial) & goCues < times.endTime(trial));

            %get the codes and times for the go cues
            goCue = nan(1,numTgt);
            if isempty(idxGo)
                tgtsAttempted = 0;
            else
                tgtsAttempted = length(idxGo);
            end
            if tgtsAttempted>0
                goCue(1:tgtsAttempted)=goCues(idxGo);
            end
            % extract first go cue
            targStart = goCue(1);

            %identify trials with corrupt end codes that might end up with extra
            %targets
            if length(idxGo) > numTgt
                warning('trt_trial_table: Inconsistent number of targets @ t = %.3f, skipping trial:%d',times.startTime(trial),trial);
                corruptDB=1;
                continue;
            end
            %find target centers
            ctr=bytes2float(cds.databursts.db(dbidx,hdrSize+1:end));
            % Offsets, target size
            xOffset = bytes2float(cds.databursts.db(dbidx,10:13));
            yOffset = bytes2float(cds.databursts.db(dbidx,14:17));
            tgtSize = bytes2float(cds.databursts.db(dbidx,18:21));
            wsnum = bytes2float(cds.databursts.db(dbidx,22:25));

            % Build arrays
            ctHoldList(trial,:)=        ctHold;             % start time of center hold
            targStartList(trial,:)=     targStart;          % time of first target onset
            goCueList(trial,:)=         goCue;              % time stamps of go_cue(s)
            numTgts(trial)=             numTgt;             % max number of targets
            numAttempted(trial,:)=      tgtsAttempted;      % ?
            xOffsets(trial)=            xOffset;            % x offset
            yOffsets(trial)=            yOffset;            % y offset
            tgtSizes(trial)=            tgtSize;            % target size
            wsnums(trial)=              wsnum;              % workspace number
            tgtCtrs(trial,:)=           ctr;                %center positions of the targets
        end

        trials=table(ctHoldList,targStartList,goCueList,numTgts,numAttempted,xOffsets,yOffsets,tgtSizes,wsnums,tgtCtrs,...
                    'VariableNames',{'ctHoldTime','targetStartTime','goCueTime','numTgt','numAttempted','xOffset','yOffset','tgtSize','spaceNum','tgtCtr'});
        trials.Properties.VariableUnits={'s','s','s','int','int','cm','cm','cm','int','cm'};
        trials.Properties.VariableDescriptions={'time of center hold start','first target go cue time','go cue time','number of targets','number of targets attempted','x offset','y offset','target size','workspace number','target center position'};
    else
        error('rw_trial_table_hdr:BadDataburstVersion',['Trial table parsing not implemented for databursts with version: ', num2str(db_version)])
    end
    if corruptDB==1
        cds.addProblem('There are corrupt databursts with more targets than expected. These have been skipped, but this frequently relates to errors in trial table parsting with the RW task')
    end
    trials=[times,trials];
    trials.Properties.Description='Trial table for the RW task';
    %cds.setField('trials',trials)
    set(cds,'trials',trials)
    evntData=loggingListenerEventData('getRWTaskTable',[]);
    notify(cds,'ranOperation',evntData)
end
