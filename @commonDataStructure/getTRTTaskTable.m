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
    goCodes=  cds.words.word(bitand(hex2dec('f0'), cds.words.word) == wordGo)-wordGo;

    %check DB version number and run appropriate parsing code
    % DB version 0 has 17 bytes before target position
    % DB version 1 has 21 bytes before target position
    db_version=cds.databursts.db(1,2);

    if db_version==0
        hdrSize=17;
        numTgt = (cds.databursts.db(1)-hdrSize)/8;

        goCueList=      nan(numTrials,numTgt);
        goCodeList=     nan(numTrials,numTgt);
        numTgts=        numTgt*ones(numTrials,1);
        numAttempted=   nan(numTrials,1);
        xOffsets=       nan(numTrials,1); 
        yOffsets=       nan(numTrials,1);
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
                warning('rw_trial_table: Inconsistent number of targets @ t = %.3f, skipping trial:%d',times.startTime(trial),trial);
                corruptDB=1;
                continue;
            end

            % Go cues
            idxGo = find(goCues > times.startTime(trial) & goCues < times.endTime(trial));

            %get the codes and times for the go cues
            goCue = nan(1,numTgt);
            goCode= nan(1,numTgt);
            if isempty(idxGo)
                tgtsAttempted = 0;
            else
                tgtsAttempted = length(idxGo);
            end
            if tgtsAttempted>0
                goCue(1:tgtsAttempted)=goCues(idxGo);
                goCode(1:tgtsAttempted)= goCodes(idxGo);
            end

            %identify trials with corrupt end codes that might end up with extra
            %targets
            if length(idxGo) > numTgt
                warning('rw_trial_table: Inconsistent number of targets @ t = %.3f, skipping trial:%d',times.startTime(trial),trial);
                corruptDB=1;
                continue;
            end
            %find target centers
            ctr=bytes2float(cds.databursts.db(dbidx,hdrSize+1:end));
            % Offsets, target size
            xOffset = bytes2float(cds.databursts.db(dbidx,10:13));
            yOffset = bytes2float(cds.databursts.db(dbidx,14:17));

            % Build arrays
            goCueList(trial,:)=         goCue;              % time stamps of go_cue(s)
            goCodeList(trial,:)=        goCode;             % ?
            numTgts(trial)=             numTgt;             % max number of targets
            numAttempted(trial,:)=      tgtsAttempted;      % ?
            xOffsets(trial)=            xOffset;            % x offset
            yOffsets(trial)=            yOffset;            % y offset
            tgtSizes(trial)=            tgtSize;            % target size
            tgtCtrs(trial,:)=           ctr;                %center positions of the targets
        end

        trials=table(goCueList,goCodeList,numTgts,numAttempted,xOffsets,yOffsets,tgtCtrs,...
                    'VariableNames',{'goCueTime','tgtID','numTgt','numAttempted','xOffset','yOffset','tgtCtr'});
        trials.Properties.VariableUnits={'s','int','int','int','cm','cm','cm','cm'};
        trials.Properties.VariableDescriptions={'go cue time','code of the go cue','number of targets','number of targets attempted','x offset','y offset','target center position'};

    elseif dbversion==1
        hdrSize=21;
        numTgt = (cds.databursts.db(1)-hdrSize)/8;

        goCueList=      nan(numTrials,numTgt);
        goCodeList=     nan(numTrials,numTgt);
        numTgts=        numTgt*ones(numTrials,1);
        numAttempted=   nan(numTrials,1);
        xOffsets=       nan(numTrials,1); 
        yOffsets=       nan(numTrials,1);
        tgtSizes=       nan(numTrials,1); % new in dbversion 1
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
                warning('rw_trial_table: Inconsistent number of targets @ t = %.3f, skipping trial:%d',times.startTime(trial),trial);
                corruptDB=1;
                continue;
            end

            % Go cues
            idxGo = find(goCues > times.startTime(trial) & goCues < times.endTime(trial));

            %get the codes and times for the go cues
            goCue = nan(1,numTgt);
            goCode= nan(1,numTgt);
            if isempty(idxGo)
                tgtsAttempted = 0;
            else
                tgtsAttempted = length(idxGo);
            end
            if tgtsAttempted>0
                goCue(1:tgtsAttempted)=goCues(idxGo);
                goCode(1:tgtsAttempted)= goCodes(idxGo);
            end

            %identify trials with corrupt end codes that might end up with extra
            %targets
            if length(idxGo) > numTgt
                warning('rw_trial_table: Inconsistent number of targets @ t = %.3f, skipping trial:%d',times.startTime(trial),trial);
                corruptDB=1;
                continue;
            end
            %find target centers
            ctr=bytes2float(cds.databursts.db(dbidx,hdrSize+1:end));
            % Offsets, target size
            xOffset = bytes2float(cds.databursts.db(dbidx,10:13));
            yOffset = bytes2float(cds.databursts.db(dbidx,14:17));
            tgtSize = bytes2float(cds.databursts.db(dbidx,18:21));

            % Build arrays
            goCueList(trial,:)=         goCue;              % time stamps of go_cue(s)
            goCodeList(trial,:)=        goCode;             % ?
            numTgts(trial)=             numTgt;             % max number of targets
            numAttempted(trial,:)=      tgtsAttempted;      % ?
            xOffsets(trial)=            xOffset;            % x offset
            yOffsets(trial)=            yOffset;            % y offset
            tgtSizes(trial)=            tgtSize;            % target size
            tgtCtrs(trial,:)=           ctr;                %center positions of the targets
        end

        trials=table(goCueList,goCodeList,numTgts,numAttempted,xOffsets,yOffsets,tgtSizes,tgtCtrs,...
                    'VariableNames',{'goCueTime','tgtID','numTgt','numAttempted','xOffset','yOffset','tgtSize','tgtCtr'});
        trials.Properties.VariableUnits={'s','int','int','int','cm','cm','cm','cm'};
        trials.Properties.VariableDescriptions={'go cue time','code of the go cue','number of targets','number of targets attempted','x offset','y offset','target size','target center position'};

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
