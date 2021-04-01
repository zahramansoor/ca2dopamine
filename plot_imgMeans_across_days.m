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
        
%ZD's addition
        %plot mean of plane 1, plane 2, etc. for same animal across
        %multiple days
        
close all
clear all
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
%days = ["d2", "d3", "d4", "d5", "d6", "d8", "d9"]; %e156
%days = ["d2", "d3", "d5", "d6", "d7", "d8", "d9"]; %e157
days = ["d2", "d3", "d4", "d5", "d6", "d7", "d8", "d9"]; %e158
src = 'Z:\analysis\plots'; %save location for plots
%2x2 tiled image
figure('DefaultAxesFontSize',10); hold on;
for pln = 1:4 %4 = most superficial layer, 1 = deepest layer
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
        rew_idx=find(R); %get indexes of all rewards
        rew_idx_diff = diff(rew_idx); %difference in reward index from last
        short = rew_idx_diff<num_rew_win_frames; %logical for rewards that happen less than x frames from last reward. 0 = single rew.
        %------------------gerardo's edits------------------
        singlerewidx = [0 diff(short)'];
        if(singlerewidx(2) == -1)
            singelidx(1) = 1;
        end
        singlerewidx = find(singlerewidx == 0);
        single_idx = rew_idx(singlerewidx);
        %------------------gerardo's edits------------------
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
            std1 = mean(norm_double_traces,2) + std(norm_double_traces')'; %upper bound
            std2 = mean(norm_double_traces,2) - std(norm_double_traces')'; %lower bound
            
            subplot(2,2,pln)
            hold all
            title(sprintf('mean +/- SD of double rewards for plane %d', pln));
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
    end
end
%formatting to keep track of days plotted w/o legend
txt = sprintf('recording days: %s', strjoin(days, ', '));
annotation('textbox',[.9 .5 .1 .2],'String',txt,'EdgeColor','none')
currfile = strcat(src, '\', anpath(4:7), '_double_rew_across_days.fig');
savefig(currfile)
hold off

%one fig for single rewards
figure('DefaultAxesFontSize',10); hold on;
for pln = 1:4 %4 = most superficial layer, 1 = deepest layer
    for daynum = 1:length(days) %iterate through days
        clear single_idx single_lick_idx %to prevent index errors
        daypath = fullfile(anpath, days{daynum});
        string = sprintf('*mean_plane%d.mat', pln);
        matfile = dir(fullfile(daypath, string)); %find mean plane 1 for analysis
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
        rew_idx=find(R); %get indexes of all rewards
        rew_idx_diff = diff(rew_idx); %difference in reward index from last
        short = rew_idx_diff<num_rew_win_frames; %logical for rewards that happen less than x frames from last reward. 0 = single rew.

        %single rewards
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
        
        std1 = mean(norm_single_traces,2) + std(norm_single_traces')'; %upper bound
        std2 = mean(norm_single_traces,2) - std(norm_single_traces')'; %lower bound
        
        subplot(2,2,pln)
        hold all
        title(sprintf('mean +/- SD of single rewards for plane %d', pln));
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
    end
end
%formatting to keep track of days plotted w/o legend
txt = sprintf('recording days: %s', strjoin(days, ', '));
annotation('textbox',[.9 .5 .1 .2],'String',txt,'EdgeColor','none')
currfile = strcat(src, '\', anpath(4:7), '_single_rew_across_days.fig');
savefig(currfile)
hold off

%one fig for non reward licks
figure('DefaultAxesFontSize',10); hold on;
for pln = 1:4 %4 = most superficial layer, 1 = deepest layer
    for daynum = 1:length(days) %iterate through days
        daypath = fullfile(anpath, days{daynum});
        string = sprintf('*mean_plane%d.mat', pln);
        matfile = dir(fullfile(daypath, string)); %find mean plane 1 for analysis
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
        rew_idx=find(R); %get indexes of all rewards
        rew_idx_diff = diff(rew_idx); %difference in reward index from last
        short = rew_idx_diff<num_rew_win_frames; %logical for rewards that happen less than x frames from last reward. 0 = single rew.

        %non-rewarded licks
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
        
        std1 = mean(norm_nr_traces,2) + std(norm_nr_traces')'; %upper bound
        std2 = mean(norm_nr_traces,2) - std(norm_nr_traces')'; %lower bound
        
        subplot(2,2,pln)
        hold all
        title(sprintf('mean +/- SD of non-rewarded licks for plane %d', pln));
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
    end
end
%formatting to keep track of days plotted w/o legend
txt = sprintf('recording days: %s', strjoin(days, ', '));
annotation('textbox',[.9 .5 .1 .2],'String',txt,'EdgeColor','none')
currfile = strcat(src, '\', anpath(4:7), '_no_rew_across_days.fig');
savefig(currfile)  
hold off

%%
%---------------------------------------------------------------NOT EDITED BEYOND THIS POINT------------------------------------------------------------------------  
%---------------------------------------------------------------NOT EDITED BEYOND THIS POINT------------------------------------------------------------------------  
%---------------------------------------------------------------NOT EDITED BEYOND THIS POINT------------------------------------------------------------------------  
%---------------------------------------------------------------NOT EDITED BEYOND THIS POINT------------------------------------------------------------------------  
%---------------------------------------------------------------NOT EDITED BEYOND THIS POINT------------------------------------------------------------------------  
%---------------------------------------------------------------NOT EDITED BEYOND THIS POINT------------------------------------------------------------------------  
%---------------------------------------------------------------NOT EDITED BEYOND THIS POINT------------------------------------------------------------------------  
%---------------------------------------------------------------NOT EDITED BEYOND THIS POINT------------------------------------------------------------------------  
%---------------------------------------------------------------NOT EDITED BEYOND THIS POINT------------------------------------------------------------------------  
%---------------------------------------------------------------NOT EDITED BEYOND THIS POINT------------------------------------------------------------------------  
%---------------------------------------------------------------NOT EDITED BEYOND THIS POINT------------------------------------------------------------------------  
%---------------------------------------------------------------NOT EDITED BEYOND THIS POINT------------------------------------------------------------------------  
%---------------------------------------------------------------NOT EDITED BEYOND THIS POINT------------------------------------------------------------------------  
%---------------------------------------------------------------NOT EDITED BEYOND THIS POINT------------------------------------------------------------------------  
%---------------------------------------------------------------NOT EDITED BEYOND THIS POINT------------------------------------------------------------------------  
%---------------------------------------------------------------NOT EDITED BEYOND THIS POINT------------------------------------------------------------------------  
%---------------------------------------------------------------NOT EDITED BEYOND THIS POINT------------------------------------------------------------------------  
%---------------------------------------------------------------NOT EDITED BEYOND THIS POINT------------------------------------------------------------------------  
%---------------------------------------------------------------NOT EDITED BEYOND THIS POINT------------------------------------------------------------------------  
%---------------------------------------------------------------NOT EDITED BEYOND THIS POINT------------------------------------------------------------------------  

save(tifffilename,'norm_nr_traces','nr_traces','nr_lick_idx','all_rew_lick','nr_lick','-append');

 if exist('double_rew', 'var')
    figure;
    hold on; 
    title('Smoothed Mean Double Rewards, Singles, non-Rewarded Licks');
    xlabel('seconds from first reward lick')
    ylabel('dF/F')
%     plot(frame_time*(-pre_win_frames):frame_time:frame_time*post_win_frames,norm_single_traces,'Color',[.8 .8 .8]);
    plot(frame_time*(-pre_win_frames):frame_time:frame_time*post_win_frames,smoothdata(mean(norm_double_traces,2),'gaussian',gauss_win/2));
    plot(frame_time*(-pre_win_frames):frame_time:frame_time*post_win_frames,smoothdata(mean(norm_single_traces,2),'gaussian',gauss_win/2)); 
    plot(frame_time*(-pre_win_frames):frame_time:frame_time*post_win_frames,smoothdata(nanmean(norm_nr_traces,2),'gaussian',gauss_win/2)); 

%     legend(['n = ',num2str(size(norm_single_traces,2))])
    legend('Double Reward','Single Reward','Unrewarded Licks')

    currfile=strcat(stripped_tifffilename,'_doubles_singles_NR.fig');
    savefig(currfile)
    
 end

 %solenoid 2   
 if exist('sol2','var')
     sol2_binned = accumarray(bin_indx(:),sol2,[],@max);
     Sol2 = bwlabel(sol2_binned>sol2_thresh);
     sol2_rew = Sol2>0;
     sol2_idx=find(sol2_rew);
     sol2_idx_diff=diff(sol2_idx);
     sol2_short=sol2_idx_diff<num_rew_win_frames;
     %can do doubles sol2 here.  adapt doubles code from above
 
 
    if any(sol2_rew)

     figure,hold on;plot((supraLick*.01)+1); plot(((rew_binned*2)+1));plot(smoothdata(norm_base_mean,'gaussian',gauss_win));plot(smoothdata(((roe_binned/100)+1),'gaussian',gauss_win));
     plot((sol2_rew*.012)+1);
     title(['Smoothed fluorescence, licks, rewards, ROE, and Solenoid2. win= ' num2str(gauss_win)]);
     currfile=strcat(stripped_tifffilename,'_Smoothed_lick_rew_ROE_fl_sol2.fig');
     savefig(currfile)


     %single solenoid2
      sol2_multi_rew_expand=bwlabel(sol2_short);%single rewards are 0
     for i=1:length(sol2_multi_rew_expand)
           sol2_multi_rew_expand(find(sol2_multi_rew_expand==i,1,'last')+1)=i;%need to expand index of multi reward by 1 to properly match rew_ind
     end

     if length(sol2_multi_rew_expand) < length(sol2_idx)
         sol2_multi_rew_expand(end+1)=0;%need to add extra on end to match index. Above for loop does this if last rew is multi reward. this does for single last.
     end

     single_sol2=find(sol2_multi_rew_expand==0);

      for i=1:length(single_sol2)

              %single_idx(i)=rew_idx(i); %orig but doesn't eliminate doubles
             single_sol2_idx(i)=sol2_idx(single_sol2(i));
          if single_sol2_idx(i)+rew_lick_win_frames < length(sol2_rew)%if window to search for lick after rew is past length of supraLick, doesn't update single_lick_idx, but single_idx is
        %      single_lick_idx(i)= (find(supraLick(single_idx(i):single_idx(i)+num_rew_win_frames),1,'first'))+single_idx(i)-1;
              single_sol2_lick_idx(i)= (find(sol2_rew(single_sol2_idx(i):single_sol2_idx(i)+rew_lick_win_frames),1,'first'))+single_sol2_idx(i)-1;
              %looks for first lick after rew with window =exclusion_win_frames
              %however first lick can be much further in naive animals
          end
      end
    %     single_sol2_lick_idx=single_sol2_lick_idx>0;    %if can't find lick in window, remove from list of lick idx
          if single_sol2_lick_idx(1) - pre_win_frames <0%remove events too early
              single_sol2_lick_idx(1)=[];
          end
          if single_sol2_lick_idx(end) + post_win_frames > length(base_mean)%remove events too late
              single_sol2_lick_idx(end)=[];
          end
     single_sol2_traces=zeros(pre_win_frames+post_win_frames+1,length(single_sol2_lick_idx));
      for i=1:length(single_sol2_lick_idx)
         single_sol2_traces(:,i)=base_mean(single_sol2_lick_idx(i)-pre_win_frames:single_sol2_lick_idx(i)+post_win_frames)';%lick at pre_win_frames+1
      end
        norm_single_sol2_traces=single_sol2_traces./mean(single_sol2_traces(1:pre_win_frames,:));

        figure;
        hold on; 
        title('Single Solenoid2');
        xlabel('seconds from first reward lick')
        ylabel('dF/F')
    %     plot(frame_time*(-pre_win_frames):frame_time:frame_time*post_win_frames,norm_single_traces,'Color',[.8 .8 .8]);
        plot(frame_time*(-pre_win_frames):frame_time:frame_time*post_win_frames,norm_single_sol2_traces);
        plot(frame_time*(-pre_win_frames):frame_time:frame_time*post_win_frames,mean(norm_single_sol2_traces,2),'k','LineWidth',2); 
    %     legend(['n = ',num2str(size(norm_single_traces,2))])
        legend()

        currfile=strcat(stripped_tifffilename,'_single_sol2.fig');
        savefig(currfile)

        save(tifffilename,'sol2_binned','sol2_rew','sol2_idx','sol2_short','sol2_multi_rew_expand','single_sol2','single_sol2_idx','single_sol2_lick_idx','single_sol2_traces',...
            'norm_single_sol2_traces','sol2_thresh','-append')


        figure;
        hold on; 
        title('Smoothed Mean Double Rewards, Singles, non-Rewarded Licks, and Fake Rewards');
        xlabel('seconds from first reward lick')
        ylabel('dF/F')
    %     plot(frame_time*(-pre_win_frames):frame_time:frame_time*post_win_frames,norm_single_traces,'Color',[.8 .8 .8]);
        plot(frame_time*(-pre_win_frames):frame_time:frame_time*post_win_frames,smoothdata(mean(norm_double_traces,2),'gaussian',gauss_win/2));
        plot(frame_time*(-pre_win_frames):frame_time:frame_time*post_win_frames,smoothdata(mean(norm_single_traces,2),'gaussian',gauss_win/2)); 
        plot(frame_time*(-pre_win_frames):frame_time:frame_time*post_win_frames,smoothdata(nanmean(norm_nr_traces,2),'gaussian',gauss_win/2)); 
        plot(frame_time*(-pre_win_frames):frame_time:frame_time*post_win_frames,smoothdata(mean(norm_single_sol2_traces,2),'gaussian',gauss_win/2)); 
    %     legend(['n = ',num2str(size(norm_single_traces,2))])
        legend('Double Reward','Single Reward','Unrewarded Licks','Fake Reward')

        currfile=strcat(stripped_tifffilename,'doubles_singles_NR_fakes.fig');
        savefig(currfile)

    end

 end

 
 

 