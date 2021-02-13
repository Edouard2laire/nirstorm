function varargout = process_nst_detect_motion( varargin )

% @=============================================================================
% This function is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2016 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Edouard Delaire, 2020; Thomas Vincent, 2015-2019


eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Detect motion artifact';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = {'NIRS', 'Pre-process'};
    sProcess.Index       = 1300; %0: not shown, >0: defines place in the list of processes
    sProcess.Description = 'http://neuroimage.usc.edu/brainstorm/Tutorials/NIRSFingerTapping#Bad_channel_tagging';
    sProcess.isSeparator = 0; 
    
    sProcess.options.windows_length.Comment = 'Windows length';
    sProcess.options.windows_length.Type    = 'value';
    sProcess.options.windows_length.Value   = {10,'s',0};

    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data', 'raw'};
    % Definition of the outputs of this process
    sProcess.OutputTypes = {'data', 'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    
   

end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
        
    % Load channel file
    ChanneMat = in_bst_channel(sInputs(1).ChannelFile);
    
    % Load recordings
    if strcmp(sInputs.FileType, 'data')     % Imported data structure
        sData = in_bst_data(sInputs(1).FileName);
    elseif strcmp(sInputs.FileType, 'raw')  % Continuous data file       
        sData = in_bst(sInputs(1).FileName, [], 1, 1, 'no');
    end
   
    sProcess.options.windows_length=sProcess.options.windows_length.Value{1}*10;
    [TOI,motion] = Compute(sData, ChanneMat, sProcess.options);
    
    nirs_ichans = (sData.ChannelFlag==1)'  & strcmpi({ChanneMat.Channel.Type}, 'NIRS');
    sEvent = db_template('event');
    sEvent.times=motion;
    sEvent.epochs=ones(1,size(motion,2));
    sEvent.label='motion';
    sData.Events=[sData.Events sEvent];
    
    
    bst_save(file_fullpath(sInputs(1).FileName), sData, 'v7');
    

    sDataOut                    = sData;
    sDataOut.F(nirs_ichans,:)   = TOI';
    sDataOut.Comment            = [sInputs.Comment '| TOI'];
    sDataOut                    = bst_history('add', sDataOut, 'process_nst_detect_motion','detect motion');

    % Generate a new file name in the same folder
    sStudy = bst_get('Study', sInputs.iStudy);
    OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_toi');
    sDataOut.FileName = file_short(OutputFile);
    bst_save(OutputFile, sDataOut, 'v7');
    % Register in database
    db_add_data(sInputs.iStudy, OutputFile, sDataOut);
    OutputFiles{1} = OutputFile;

    
end


%% ===== Compute =====
function [TOI,motion]= Compute(sData, channel_def, options)

    % Todo: 
    % - Metric: check the difference using different metric: std, sci, power, 
    % - Distribution fit: see the importance of the bandwith parameter
    % - Votation: determine a better threeshold value. (naive bayes
    % classifier?)
    
    nirs_flags = (sData.ChannelFlag==1)'  & strcmpi({channel_def.Channel.Type}, 'NIRS');
    channel_def = channel_def.Channel(nirs_flags);
    windows_length = options.windows_length;
    
    signal=sData.F';
    n_sample = length(sData.Time);
    nirs_signal=signal(:,nirs_flags);
    n_channel = size(nirs_signal,2);
    
    mov_mean= zeros(n_sample,n_channel); 
    mov_std= zeros(n_sample,n_channel); 
    
    
    % Note: the following code is equivalent to : 
    % mov_mean = movmean(nirs_signal,windows_length+1,1);
    % mov_std  = movstd(nirs_signal,windows_length+1,1);
    
    for i = 1:n_sample
       windows_start =  max( 1,i - windows_length/2);
       windows_end   =  min( length(sData.Time), i + windows_length/2);
       
       nirs_windows  = nirs_signal(windows_start:windows_end,:);
       mov_mean(i,:) = mean(nirs_windows);
       mov_std(i,:)  = std(nirs_windows);
        
    end    
    
    CV = mov_std./mov_mean.*100;
    
    threshold=zeros(1,n_channel);
    TOI= zeros(size(nirs_signal)); % Time of interest.
    
    figure(2); clf;
    subplot(ceil(sqrt(n_channel)),ceil(sqrt(n_channel)),1);
    warning('')
    
    for i_chan=1:n_channel               
        [f,xi] = ksdensity(CV(:,i_chan),'Bandwidth',0.01);
        [est_icdft,xi2,bw] = ksdensity(CV(:,i_chan),'Function','icdf'); %,'Bandwidth',0.01);
        [warnMsg, warnId] = lastwarn;
        disp(sprintf('channel %s - BW: %.2f',channel_def(i_chan).Name,bw)) 
        if ~isempty(warnMsg)
            disp(sprintf('icdf did not converge for channel %s (%s)',...
                         channel_def(i_chan).Name, warnMsg));
        end     
   
        if isempty(min(est_icdft(xi2 >= 0.95)))
            threshold(i_chan) = max(CV(:,i_chan)) +1;
        else
            threshold(i_chan) = min(est_icdft(xi2 >= 0.95));   
        end
        TOI(:,i_chan)= abs(CV(:,i_chan)) > threshold(i_chan);

        subplot(ceil(sqrt(n_channel)),ceil(sqrt(n_channel)),i_chan);
        plot(xi,f)
        line([threshold(i_chan) threshold(i_chan)], ylim(gca),'Color', 'r') 
        title(  channel_def(i_chan).Name)      
    end
    
    vote=sum(TOI,2); 
    % percentage_of_vote=sum(TOI,1); 
    % arround 5% for all channels as expected :)
    % figure; stem(100*percentage_of_vote/n_sample)
    
    % Threshold the vote using hysteresis principle
    threshold_high = 12;
    threshold_low = 5;
    
    is_motion=0;
    motion=[];
    for i= 1:length(vote)
       if  vote(i) > threshold_high && ~is_motion
            start_motion=sData.Time(i);
            is_motion=1;
       end
       
       if is_motion && vote(i) < threshold_low
            is_motion=0;
            end_motion = sData.Time(i);
            
            motion(:,end+1)=[start_motion ; end_motion];
       end    
    end    
    
    figure;
    subplot(2,1,1);
    hold on;    
    plot(sData.Time,nirs_signal);title('NIRS')
    y_max=max(ylim(gca));
    for i_motion=1:size(motion,2)
        rectangle('Position',[motion(1,i_motion) 0  motion(2,i_motion)-motion(1,i_motion) y_max], ...
            'FaceColor',[1 0 0 0.5],'EdgeColor', [1, 0, 0, 0])
    end    
    
    subplot(2,1,2);hold on; plot(sData.Time,vote);title('Votation');
    plot(sData.Time, threshold_high*ones(size(sData.Time)),'r')
    plot(sData.Time, threshold_low*ones(size(sData.Time)),'b')

    y_max=max(ylim(gca));
    for i_motion=1:size(motion,2)
        rectangle('Position',[motion(1,i_motion) 0   motion(2,i_motion)-motion(1,i_motion) y_max], ...
            'FaceColor',[1 0 0 0.5],'EdgeColor', [1, 0, 0, 0])
    end    
    
    
    
    [B,I] = sort(mean(nirs_signal,1));
    channel_def =channel_def(I);

    figure;imagesc(TOI(:,I)')
    
%     figure;
%     subplot(4,1,1);plot(sData.Time,nirs_signal);title('NIRS')
%     subplot(4,1,2);plot(sData.Time,mov_mean);title('Mov mean')
%     subplot(4,1,3);plot(sData.Time,mov_std);title('Mov std')
%     subplot(4,1,4);plot(sData.Time,CV);title('CV')
%     
%     figure;
%     subplot(4,1,1);plot(sData.Time,nirs_signal(:,chan_number));title('NIRS')
%     ylim([0 max(nirs_signal(:,chan_number))*1.2]);
%     subplot(4,1,2);plot(sData.Time,mov_mean(:,chan_number));title('Mov mean')
%     subplot(4,1,3);plot(sData.Time,mov_std(:,chan_number));title('Mov std')
%     subplot(4,1,4);plot(sData.Time,CV(:,chan_number));title('CV')
        
end
