clear all;
load('Z:\e157\d6\behavior\vr\E157_27_Mar_2021_time(14_32_20).mat')
figure('DefaultAxesFontSize',16)
hAxis(1) = subplot(3,1,1);
plot(VR.reward(20000:30000),'LineWidth',5)
set(gca,'xtick',[])
ylabel('reward')
xlim([0,10000])
hAxis(2) = subplot(3,1,2);
plot(VR.ROE(20000:30000)*-1)
set(gca,'xtick',[])
ylabel('locomotion')
xlim([0,10000])
ylim([0,200])
hAxis(3) = subplot(3,1,3);
plot(VR.lick(20000:30000))
ylabel('lick')
set(gca,'ytick',[])
xlabel('seconds')
xlim([0,10000])
