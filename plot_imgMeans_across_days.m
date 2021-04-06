%Ed's markup:
        %plot and analyze data from GRAB-DA1h, 
        %makes peri-rewarded lick plots of mean image intensity.
        %split into single, double, and three or more rewards. searches for first
        %lick after reward to align.
        %also plots non-rewarded licks that are outside of defined lick window
        %used after concatenate_tif_mean_EH and abffileSelectStartEnd_mean_img

        %to do
        %make sure triples are plotted and saved for older experiments
        %start and stop trigger ave
        %statistical tests of intensity change
        %fit decay
        %could also use template search for non-reward related DA increases.
        %currently deletes rewards too close to end of beginning, instead could pad
        %with NaN
        %chunk licks into near and far from reward
        %track time since last reward or last double
        %rolling average of rewards
        
%ZD's additions
        %plot mean of plane 1, plane 2, etc. for same animal across
        %multiple days
        %NOTE: does not include solenoid plots since no fake rewards were
        %given
        
close all; clear all
[tifffilename, tiffpath] = uigetfile('*.mat','pick mat file'); %get path to one of the mat files for 1 day / 1 animal
cd (tiffpath); %set path
%get animal's folder path; i.e. up one folder
idcs = strfind(tiffpath,'\');
anpath = tiffpath(1:idcs(end)-4);
numplanes = 4; %ZD defined this since before it was loaded w the mat file

%------------------ed's params - unchanged------------------
gauss_win=12;
frame_rate=31.25/numplanes;
lickThresh=-0.07; %ZD changed bc code was crashing otherwise...
rew_thresh=0.001;
sol2_thresh=1.5;
num_rew_win_sec=5;%window in seconds for looking for multiple rewards
rew_lick_win=20;%window in seconds to search for lick after rew. could be long in naive animals but likely short latency in trained
pre_win=5;%pre window in s for rewarded lick average
post_win=20;%post window in s for rewarded lick average
exclusion_win=20;%exclusion window pre and post rew lick to look for non-rewarded licks
%------------------ed's params - unchanged------------------

%need to load 1 plane at a time per day to collect this info
%only care about these days (MAY NEED TO MODIFY LIST IF MORE DAYS ARE
%NEEDED)
%the way i'm doing this is unnecessary, figure out how to put in one dict???
%one fig for double rewards
days = ["d2", "d3", "d4", "d5", "d6", "d8", "d9", "d10", "d11", "d12"]; %e156, need to skip d7
%bc no clampex file
%days = ["d2", "d3", "d5", "d6", "d7", "d8", "d9", "d10", "d11", "d12"]; %e157, need to skip d4
%unidirectional day
%days = ["d2", "d3", "d4", "d5", "d6", "d7", "d8", "d9", "d10", "d11"]; %e158
src = 'Z:\analysis\plots'; %save location for plots
%open 3 figures
%2x2 tiled image
fig1 = figure('DefaultAxesFontSize',10); hold on; %double rew
fig2 = figure('DefaultAxesFontSize',10); hold on; %single rew
fig3 = figure('DefaultAxesFontSize',10); hold on; %no rew
fig4 = figure('DefaultAxesFontSize',10); hold on; %double, single, no rew overlay
fig5 = figure('DefaultAxesFontSize',10); hold on; %fluo during 2 min delay period
for pln = 1:numplanes %4 = most superficial layer, 1 = deepest layer
    for daynum = 1:length(days) %iterate through days
        daypath = fullfile(anpath, days{daynum});
        string = sprintf('*mean_plane%d.mat', pln);
        matfile = dir(fullfile(daypath, string)); %find mean plane #num for analysis
        matfl = fullfile(daypath, matfile.name); %join path to mat file for the plane
        disp(matfl)
        disp(' ')
        load(matfl);
        %------------------ed's calculated values - unchanged------------------
        frame_time = 1/frame_rate;
        num_rew_win_frames = round(num_rew_win_sec/frame_time);%window in frames
        rew_lick_win_frames = round(rew_lick_win/frame_time);%window in frames
        post_win_frames = round(post_win/frame_time);
        pre_win_frames = round(pre_win/frame_time);
        exclusion_win_frames = round(exclusion_win/frame_time);
        [B,~,bin_indx] = histcounts(1:numframes,length(base_mean));
        rew_binned = accumarray(bin_indx(:),rewards,[],@mean);
        mean_base_mean = mean(base_mean);
        
        norm_base_mean = base_mean/mean_base_mean;
        lick_binned = accumarray(bin_indx(:),lick,[],@min);
        roe_binned = accumarray(bin_indx(:),ROE,[],@max);
        L = bwlabel(lick_binned<lickThresh);
        supraLick = L > 0;

        mean_lick = mean(lick_binned);
        %------------------ed's calculated values - unchanged------------------

        %ed's vars for reward location + analysis
        R = bwlabel(rew_binned>rew_thresh); %label rewards, ascending
        rew_idx = find(R); %get indexes of all rewards
        rew_idx_diff = diff(rew_idx); %difference in reward index from last
        short = rew_idx_diff<num_rew_win_frames; %logical for rewards that happen less than x frames from last reward. 0 = single rew, 1 = double rew
        %------------------gerardo's edits------------------
        singlerewidx = [0 diff(short)'];
        if(singlerewidx(2) == -1)
            singlerewidx(1) = 1;
        end
        singlerewidx = find(singlerewidx == 0);
        single_idx = rew_idx(singlerewidx);
        %------------------gerardo's edits------------------
        %DOUBLE REWARDS       
        set(0,'CurrentFigure',fig1)
        %if there are any multi rewards
        if any(short)
            multi_reward_num = bwlabel(short);%label all multiple rewards ascending. doubles have single number, triples two consecutive, etc.
            % double_rew=find(ysize==0);
            for i=1:max(multi_reward_num)  %find number of consecutive rewards < window. 
                ysize(i) = find(multi_reward_num==i,1,'last')-find(multi_reward_num==i,1,'first');
                %ysize length is number of multirewards,
            end

            %double rewards. must be doubles is any(short), assuming not just tri.
            double_rew = find(ysize==0);%double events have ysize=0
            for i=1:length(double_rew)
                double_idx(i) = rew_idx(find(multi_reward_num==double_rew(i)));
                if double_idx(i)+rew_lick_win_frames < length(supraLick)%if window to search for lick after rew is past length of supraLick, doesn't update single_lick_idx, but single_idx is
                    double_lick_idx(i) = (find(supraLick(double_idx(i):double_idx(i)+rew_lick_win_frames),1,'first'))+double_idx(i)-1;%finding closest lick after rew
                end
            end

            %remove trials too close to beginning or end
            if double_lick_idx(1) - pre_win_frames <0%remove events too early
                double_lick_idx(1)=[];
            end
            if double_lick_idx(end) + post_win_frames > length(base_mean)%remove events too late
                double_lick_idx(end)=[];
            end
            double_traces=zeros(pre_win_frames+post_win_frames+1,length(double_lick_idx));
            for i = 1:length(double_lick_idx)
                double_traces(:,i)=base_mean(double_lick_idx(i)-pre_win_frames:double_lick_idx(i)+post_win_frames);%lick at pre_win_frames+1
            end
            
            norm_double_traces=double_traces./mean(double_traces(1:pre_win_frames,:));
            std1 = prctile(norm_double_traces', 25)'; %25 p 
            std2 = prctile(norm_double_traces', 75)'; %75 p
            
            subplot(2,2,pln)
            hold all
            title(sprintf('mean, 25 & 75%ile of double rewards for plane %d', pln));
            xlabel('seconds from reward lick')
            ylabel('dF/F')
            x = frame_time*(-pre_win_frames):frame_time:frame_time*post_win_frames;
            line = plot(x, mean(norm_double_traces,2),...
                'LineWidth',2,'DisplayName',days{daynum}); %autocolor
            x2 = [x, fliplr(x)];
            inbtw = [std1', fliplr(std2')]; %shading std
            h = fill(x2, inbtw, get(line, 'Color'), 'LineStyle', 'none');%);
            set(h,'facealpha',.2);
            clear ysize %to prevent index errors
            hold off
        end
                
        %SINGLE REWARDS
        set(0,'CurrentFigure',fig2)
        multi_rew_expand=bwlabel(short);%single rewards are 0
        for i=1:length(multi_rew_expand)
            multi_rew_expand(find(multi_rew_expand==i,1,'last')+1)=i;%need to expand index of multi reward by 1 to properly match rew_ind
        end

        if length(multi_rew_expand) < length(rew_idx)
            multi_rew_expand(end+1)=0;%need to add extra on end to match index. Above for loop does this if last rew is multi reward. this does for single last.
        end
        single_rew=find(multi_rew_expand==0);
        for i=1:length(single_rew)
            single_idx(i) = rew_idx(single_rew(i));
            if single_idx(i)+rew_lick_win_frames < length(supraLick) %if window to search for lick after rew is past length of supraLick, doesn't update single_lick_idx, but single_idx is
                try %ZD added for GRAB-D1Ah dopamine mice that sometimes don't lick on the first try
                    single_lick_idx(i)= (find(supraLick(single_idx(i):single_idx(i)+rew_lick_win_frames),1,'first'))+single_idx(i)-1;
                catch 
                    warning('lick not detected around reward, likely the animal missed it. skipping index...')
                end
                %looks for first lick after rew with window =exclusion_win_frames
                %however first lick can be much further in naive animals
            end
        end

        if single_lick_idx(1) - pre_win_frames <0%remove events too early
            single_lick_idx(1)=[];
        end
        if single_lick_idx(end) + post_win_frames > length(base_mean)%remove events too late
            single_lick_idx(end)=[];
        end
        single_lick_idx = nonzeros(single_lick_idx)'; %removes 0's assigned when lick wasn't detected around reward

        single_traces = zeros(pre_win_frames+post_win_frames+1,length(single_lick_idx));
        for i=1:length(single_lick_idx)
            single_traces(:,i) = base_mean(single_lick_idx(i)-pre_win_frames:single_lick_idx(i)+post_win_frames)';%lick at pre_win_frames+1
        end
        norm_single_traces=single_traces./mean(single_traces(1:pre_win_frames,:));
        
        std1 = prctile(norm_single_traces', 25)'; %25 p 
        std2 = prctile(norm_single_traces', 75)'; %75 p 
        
        subplot(2,2,pln)
        hold all
        title(sprintf('mean, 25 & 75%ile of single rewards for plane %d', pln));
        xlabel('seconds from first reward lick')
        ylabel('dF/F')
        x = frame_time*(-pre_win_frames):frame_time:frame_time*post_win_frames;
        line = plot(x, mean(norm_single_traces,2),...
            'LineWidth',2,'DisplayName',days{daynum}); %autocolor
        x2 = [x, fliplr(x)];
        inbtw = [std1', fliplr(std2')]; %shading std
        h = fill(x2, inbtw, get(line, 'Color'), 'LineStyle', 'none');%);
        set(h,'facealpha',.1);
        hold off      
              
        %NON-REWARDED LICKS
        set(0,'CurrentFigure',fig3)
        all_rew_lick = single_lick_idx;
        if exist('double_lick_idx', 'var')
            all_rew_lick = [all_rew_lick double_lick_idx];%combine arrays, could also use union but should not have replicates
        end
        if exist('triple_lick_idx', 'var')
            all_rew_lick = [all_rew_lick triple_lick_idx];
        end

        nr_lick = bwlabel(supraLick);
        for i=1:length(all_rew_lick)
            nr_lick((all_rew_lick(i)-exclusion_win_frames):(all_rew_lick(i)+exclusion_win_frames),1)=0;
        end
        
        nr_lick(1:exclusion_win_frames,1)=0;%get rid of non-rewarded licks at start, otherwise will crash when grab traces
        nr_lick(end-exclusion_win_frames:end,1)=0;%same at end
        nr_lick=bwlabel(nr_lick);

        for i=1:max(nr_lick)
            nr_lick_idx(i) = (find(nr_lick==i,1,'first'));%
        end
        nr_traces=zeros(pre_win_frames+post_win_frames+1,length(single_lick_idx));
        for i=1:length(nr_lick_idx)
            nr_traces(:,i)=base_mean(nr_lick_idx(i)-pre_win_frames:nr_lick_idx(i)+post_win_frames)';%lick at pre_win_frames+1
        end
        norm_nr_traces = nr_traces./mean(nr_traces(1:pre_win_frames,:));
        
        std1 = prctile(norm_nr_traces', 25)'; %25 p 
        std2 = prctile(norm_nr_traces', 75)'; %75 p 
        
        subplot(2,2,pln)
        hold all
        title(sprintf('mean, 25 & 75%ile of non-rewarded licks for plane %d', pln));
        xlabel('seconds from non-rewarded lick')
        ylabel('dF/F')
        x = frame_time*(-pre_win_frames):frame_time:frame_time*post_win_frames;
        line = plot(x, mean(norm_nr_traces,2),...
            'LineWidth',2,'DisplayName',days{daynum}); %autocolor
        x2 = [x, fliplr(x)];
        inbtw = [std1', fliplr(std2')]; %shading std
        h = fill(x2, inbtw, get(line, 'Color'), 'LineStyle', 'none');%);
        set(h,'facealpha',.1);
        hold off 
        
        %plot double, single, no reward mean overlays
        %does not have different colors/line for the plot
        if exist('double_rew', 'var')
            set(0,'CurrentFigure',fig4)
            subplot(2,2,pln)
            hold all
            title(sprintf('smoothed mean double & single rewards & non-rewarded licks for plane %d', pln));
            xlabel('seconds from first reward lick')
            ylabel('dF/F')
            plot(frame_time*(-pre_win_frames):frame_time:frame_time*post_win_frames,...
                smoothdata(mean(norm_double_traces,2),'gaussian',gauss_win/2),'Color','r');
            plot(frame_time*(-pre_win_frames):frame_time:frame_time*post_win_frames,...
                smoothdata(mean(norm_single_traces,2),'gaussian',gauss_win/2),'Color','b'); 
            plot(frame_time*(-pre_win_frames):frame_time:frame_time*post_win_frames,...
                smoothdata(nanmean(norm_nr_traces,2),'gaussian',gauss_win/2),'Color',[0.4 0.4 0.4]);     
            hold off
        end
       
        %MEAN FLUORESCENCE DURING 1ST 2 MIN (DELAY PERIOD)
        %find indices corresponding to first 2 min
        [sz1,sz2] = size(norm_base_mean);
        indstart = sz2/(50000/31.25) * 240; %50000/31.25 gives us total imaging time in s
        fluostart = base_mean(1:ceil(indstart));
        norm_fluostart = base_mean(1:ceil(indstart))/mean(base_mean(1:ceil(indstart)));
        set(0,'CurrentFigure',fig5)
        subplot(2,2,pln)
        hold all
        title(sprintf('mean fluorescence in first 240 s, plane %d', pln));
        xlabel('frames from start of session')
        ylabel('dF/F')
        plot(norm_fluostart, 'LineWidth',0.5)
        hold off
    end
end       
%fig for double rewards
%formatting to keep track of days plotted w/o legend
hold all
txt = sprintf('recording days: %s', strjoin(days, ', '));
annotation('textbox',[.9 .5 .1 .2],'String',txt,'EdgeColor','none')
currfile = strcat(src, '\', anpath(4:7), '_double_rew_across_days.fig');
savefig(fig1, currfile)
hold off
%fig for single rewards
txt = sprintf('recording days: %s', strjoin(days, ', '));
annotation('textbox',[.9 .5 .1 .2],'String',txt,'EdgeColor','none')
currfile = strcat(src, '\', anpath(4:7), '_single_rew_across_days.fig');
savefig(fig2, currfile)
hold off
%one fig for non reward licks
%formatting to keep track of days plotted w/o legend
txt = sprintf('recording days: %s', strjoin(days, ', '));
annotation('textbox',[.9 .5 .1 .2],'String',txt,'EdgeColor','none')
currfile = strcat(src, '\', anpath(4:7), '_no_rew_across_days.fig');
savefig(fig3, currfile)  
hold off
%fig for reward overlays
txt = sprintf('recording days: %s', strjoin(days, ', '));
annotation('textbox',[.9 .5 .1 .2],'String',txt,'EdgeColor','none')
currfile = strcat(src, '\', anpath(4:7), '_double_single_no_rew_across_days.fig');
savefig(fig4, currfile)  
hold off
%fig for first 120s 
txt = sprintf('recording days: %s', strjoin(days, ', '));
annotation('textbox',[.9 .5 .1 .2],'String',txt,'EdgeColor','none')
currfile = strcat(src, '\', anpath(4:7), '_first_2min_delay_period.fig');
savefig(fig4, currfile)  
hold off
        