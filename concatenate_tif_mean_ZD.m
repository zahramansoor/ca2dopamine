%developed for GRAB-DA1h analysis, with FOV fluorescence changes
%concatenates all tifs in dir, splits into planes, baseline corrects, gets
%mean image intensity, and saves tifs of baseline corrected planes.
%also plots raw mean and baselined mean
%used after runVideosTiff_EH_new_sbx_uint16 or runVideosTiff_EH_new_sbx
%followed by abffileSelectStartEnd_mean_img and plot_imgMeans_rew_licks
%no motion correction in this path
%saves raw mean image fluorescence and baseline figs for each plane

%also works on motion stabilized output of suite2p but timebase off
%make sure no other tifs in dir

%to do

%could also bypass runVideos but need to rewrite
%could crop image prior to baseline, might help with motion artifact for
%planes that don't have fluorescence filling the field

close all
clear all

img_write=1;%1 to save baseline corrected tifs, 0 to skip

numplanes=4;
%crop coordinates. full size after runVideos (bi), x:629, y=512
%crop top to remove ETL artifact
%uncomment to crop
x1=41; % ZD changed this dim to remove stripy lines
x2=619;
y1=169;
y2=512;

base_window=200;
lenVid=2000;
[tifffilename,tiffpath]=uigetfile('*.tif','pick your tif file');
stripped_tifffilename=regexprep(tifffilename,'.tif','');  
%tic
cd (tiffpath); %set path

tifFiles = dir('*.tif'); %reads sequence all jumbled *-1, *-10, etc)
num_files = length(tifFiles);
% mydata = cell(1, numfiles);
cated_movie=[];
cated_mean=[];
chone=[];

for k=1:numplanes
    for j=1:num_files

        myfilename = (tifFiles(j).name);
        info=imfinfo(myfilename);
        numframes=length(info);
        M=info(1).Width;
        N=info(1).Height;

        chone_temp=zeros(N,M,numframes);
%         mean_temp=zeros(numframes);
        for i=1:numframes
            chone_temp(:,:,i)=imread(myfilename,i,'Info',info);
            
        end

         chone_temp=chone_temp(:,:,k:numplanes:end);   
        chone=cat(3,chone,chone_temp);
    end
        clear chone_temp
        numframes_all = size(chone,3);
        
        if exist('x1','var')
            chone=chone(max(y1,1):min(y2,N),max(x1,1):min(x2,M),:);%
        end
                
        raw_mean=zeros(1,numframes_all);
        for i=1:numframes_all
            raw_mean(i)=mean(mean(chone(:,:,i)));
        end
        figure,plot(raw_mean);
        title('Image Raw Fluorescence')
        currfile=strcat(stripped_tifffilename,'_mean_plane%d_raw.fig');
        currfile=sprintf(currfile,k);
        savefig(currfile)
        

        meanlastframes=median(mean(mean(chone(:,:,(end-base_window):end))));
        meanfirstframes=median(mean(mean(chone(:,:,1:base_window))));
        chone=chone*(meanlastframes/meanfirstframes);
        %baseline subtract whole movie
        junk=squeeze(mean(mean(chone)));%mean of each image, frame x 1 vector
        mean_all=mean(mean(mean(chone)));
        junk2=zeros(size(junk));
            for kk=1:length(junk)
                cut=junk(max(1,kk-base_window):min(numframes_all,kk+base_window));
                cutsort=sort(cut);
                a=round(length(cut)*.08);
                junk2(kk)=cutsort(a);
            end

            for i=1:numframes_all
                chone(:,:,i)=(chone(:,:,i)/junk2(i))*mean_all;
            end
        chone=uint16(chone);

        base_mean=zeros(1,numframes_all);
        for i=1:numframes_all
            base_mean(i)=mean(mean(chone(:,:,i)));
        end
        figure,plot(base_mean);
        title('Baseline Corrected Fluorescence')
        currfile=strcat(stripped_tifffilename,'_mean_plane%d_base.fig');
        currfile=sprintf(currfile,k);
        savefig(currfile)
        
        
     currfile=strcat(stripped_tifffilename,'_mean_plane%d.mat');
        currfile=sprintf(currfile,k);
        currfilename=[tiffpath currfile];
        save(currfilename,'raw_mean','base_mean','base_window','numplanes');    %need -v7.3 MAT file or variable is too big to save
        
   currfile=regexprep(currfile,'.mat','');
    if img_write==1
        for ii=1:ceil(numframes_all/lenVid) %splitting into  frame chunks. ii=1:number of files
        % ii=1;
            if ii>9
                currfile_new=strcat(currfile,'_x',num2str(ii),'.tif');
            else
                currfile_new=strcat(currfile,'_',num2str(ii),'.tif');
            end

            if ii*lenVid > size(chone,3)
                chtemp=chone(:,:,((ii-1)*lenVid)+1:end);
            else
                chtemp=chone(:,:,((ii-1)*lenVid)+1:ii*lenVid);
            end

        %         chtemp=chtemp(:,90:718,:);
            final_filename=[tiffpath currfile_new]; %files need to have path
               imageJ_savefilename=strrep(final_filename,'\','\\'); %ImageJ needs double slash
            imageJ_savefilename=['path=[' imageJ_savefilename ']'];
            % Miji;    %calls Fiji
        %     MIJ.createImage('chone_image', gray2ind(mat2gray(chtemp,[0 32767])), true);
        %     MIJ.createImage('chone_image', gray2ind(mat2gray(chtemp,double(lims)),double(round(ceil(lims(2))/2))), true);
            MIJ.createImage('chone_image', uint16(chtemp), true); %creates ImageJ file with 'name', matlab variable name

        %     MIJ.createImage('chone_image', int16(chtemp), true); %creates ImageJ file with 'name', matlab variable name
            MIJ.run('Save', imageJ_savefilename);   %saves with defined filename
            MIJ.run('Close All');
            
            cd (tiffpath); %set path

        end
    end
    
        
    chone=[];
end

MIJ.exit;
