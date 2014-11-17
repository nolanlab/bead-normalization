function normalize_folder(r)



global figpos
global num_beads

try   
    
    [figpos,questpos]=setFigurePositions;
    
    [data_dir, f,num_files]=getFileList;
    filenames=strcat(data_dir,f);
    
    timeOrder=orderFilenames(filenames);
    filenames=filenames(timeOrder);
    f=f(timeOrder);
    
    % -100 subtraction 
    sbtrct=false;
    
    %default if there is not an existing log
    logtxt=[];
    
    %make normalization folder if there isn't one
    norm_dir=[data_dir 'normed' filesep];
    if ~exist(norm_dir,'dir')
        mkdir(norm_dir)
        %open norm folder and read log
    elseif exist([norm_dir 'normalizer.log'],'file')
        oldlog=fopen([norm_dir 'normalizer.log'],'r');
        logtxt=textscan(oldlog,'%s','Delimiter','\n');
        logtxt=logtxt{1};
    end
    
    %start writing log
    fid=fopen([norm_dir 'normalizer.log'],'w');
    fprintf(fid,'%s\n',datestr(now));
    fprintf(fid,'%s\n',['Normalizing files in folder ' data_dir]);
    
    %identify beads from all files
    [bead_metals, bead_data, normthisfile]=identifyBeads;

    % beads removal option
    rmv=false;
    if isdeployed
        quest=questdlg('Do you want to remove beads?','Bead Removal','Yes','No','Yes');
        if strcmp(quest,'Yes')
            rmv=true;
        end
    elseif nargin>0
         rmv=r;   
    end
    
    %bead directory is inside normed directory
    bead_dir=[data_dir 'beads' filesep];
    if ~exist(bead_dir,'dir')
        mkdir(bead_dir)
    end
    
    %update to not include skipped files
    num_files=nnz(normthisfile);
    if num_files==0 %only empty files were in folder
        fclose(fid);
        return
    end
    filenames=filenames(normthisfile);
    f=f(normthisfile);
    bead_data=bead_data(normthisfile);
    
    allbeads=cat(1,bead_data{:});
    
    %normalize to current files or to previously saved files?
    bead_means=setNormalizationBaseline;
    
    %smoothing window
    bead_window=200;
    
    %normalize events
   [normedbeads,smoothed,normedsmoothed]=normalizeFiles;
    
    %adjust times for plotting after concatenation (don't need with new
    %plotting)
%     bead_data=concatenateTimestamps(bead_data);
%     normedbeads=concatenateTimestamps(normedbeads);
    
    plotBeforeAndAfter

    
catch err
    if isdeployed
        str=sprintf('%s\n%s\n%s%s%g',err.identifier, err.message, err.stack.name, ' line ',err.stack.line);
        h=warndlg(str ,'Error');
        waitfor(h)
        return
    else
        rethrow(err)
    end
end

    function [bead_metals, bead_data, normthisfile]=identifyBeads
        
        %% do the bead identification in all the files
        xp=[];
        yp=[];
        normthisfile=true(num_files,1);
        bead_data=cell(1,num_files);
        for i=1:num_files
            file=norm_cytof(filenames{i},sbtrct);
            [~,name,ext]=fileparts(filenames{i});
            nameExt=[name ext];
            
            %if fcs file had no events, skip it
            if isempty(file.data)
                h=warndlg(['Skipping file ' nameExt ' because it is empty!'] ,'Empty file!');
                normthisfile(i)=false;
                waitfor(h)
                drawnow
                continue
            end
            
            %if fcs file does not have monotonically increasing time column, skip
            %it
            if any(diff(file.data(:,file.time_col))<0)
                h=warndlg(['Skipping file ' nameExt ' because its time column is not monotonically increasing!'] ,'Bad concatenation!');
                normthisfile(i)=false;
                waitfor(h)
                drawnow
                continue
            end
            
            %establish which metals we are using only once
            if i==1 || ~exist('bead_metals','var') %first file could have been empty
                file=file.find_bead_masses;
                bead_metals=file.bead_metals;
                if length(bead_metals)<2
                    h=warndlg('You must choose at least two bead channels!' ,'Bead Channels Not Selected');
                    waitfor(h)
                    drawnow
                    file=file.find_bead_masses;
                    bead_metals=file.bead_metals;
                    if length(bead_metals)<2
                        fclose(fid);
                        return
                    end
                end
            else
                file=file.find_bead_channels(bead_metals);
                %skip files missing selected bead channels
                if any(file.bead_channels==0)
                    h=warndlg(['Skipping file ' nameExt ' because it is missing the following bead channels: ' bead_metals(file.bead_channels==0) '!'] ,'Missing selected bead channel!');
                    normthisfile(i)=false;
                    waitfor(h)
                    drawnow
                    continue
                end
            end
            
            
            %if no beads found, get another chance or skip norming file
            while 1
                if ~isempty(logtxt)
                    file_ind=find(strcmp(nameExt,logtxt));
                    if ~isempty(file_ind)
                        xp=cell(1,num_beads);
                        yp=cell(1,num_beads);
                        bead_ind=file_ind+2;
                        for j=1:num_beads
                            [~,bnds]=regexp(logtxt{bead_ind},'\[(.+) (.+)\]','match','tokens');
                            xp{j}([1 4 5])=asinh(1/5*str2double(bnds{1}{1}));
                            xp{j}([2 3])=asinh(1/5*str2double(bnds{1}{2}));
                            bead_ind=bead_ind+1;
                        end
                        [~,bnds]=regexp(logtxt{bead_ind},'\[(.+) (.+)\]','match','tokens');
                        for j=1:num_beads
                            yp{j}([1 2 5])=asinh(1/5*str2double(bnds{1}{1}));
                            yp{j}([3 4])=asinh(1/5*str2double(bnds{1}{2}));
                        end
                        
                        
                    end
                end
                file=file.manual_bead_identify(xp,yp);
                xp=file.xp;
                yp=file.yp;
                bead_count=sum(file.beads);
                if bead_count==0
                    quest=questdlg('No beads are in the gate.  Either try again, or skip normalizing this file.','No Beads Found!','Re-gate','Skip file','Re-gate');
                    if strcmp(quest,'Skip file')
                        normthisfile(i)=false;
                        break
                    end
                else
                    break
                end
            end
            
            if ~normthisfile(i)
                continue
            end
            
            perc_beads=100*bead_count/length(file.beads);
            
            fprintf(fid,'%s\n',nameExt);
            fprintf(fid,'\t%s\n',[num2str(bead_count) ' beads found (' num2str(perc_beads) '% of all events)']);
            for j=1:num_beads
                fprintf(fid,'\t%s%g %g%s\n',[bead_metals{j} ' bead boundary: ['],5*sinh(xp{j}(1)),5*sinh(xp{j}(2)),']');
            end
            catyp=cat(1,yp{:});
            fprintf(fid,'\t%s%g %g%s\n','DNA bead boundary: [',5*sinh(max(catyp(:,1))),5*sinh(min(catyp(:,3))),']');
            
            
            display([num2str(bead_count) ' beads found (' num2str(perc_beads) '% of all events) in file ' nameExt])
            
            file=file.assign_bead_data;
            bead_data{i}=[file.bead_data find(file.beads)]; %[time bc1 bc2 bc3 bc4 bc5 enum] UNTRANSFORMED
            
            clear file
            
        end
        
    end

    function bead_means=setNormalizationBaseline
        %set baseline to which we're normalizing
        
        quest=questdlg('Do you want to normalize to the median level of the current files, or to previously saved FCS file(s) of beads?','Normalization Baseline','Current Files','Saved Files','Current Files');
        if strcmp(quest,'Saved Files')
            [oldBead_files,oldBead_dir] = uigetfile({'*.fcs','*.FCS'},'Choose File(s) of Gated Beads','MultiSelect','on');
            if oldBead_dir==0
                h=warndlg('No bead files selected! Normalizing to the median bead level.' ,'No files selected');
                waitfor(h)
                drawnow
                bead_means=median(allbeads(:,2:(num_beads+1)));
            else
                if iscell(oldBead_files) %selected > 1 fcs file
                    
                    num_oldBeadfiles=length(oldBead_files);
                    oldBead_vals=cell(1,num_oldBeadfiles);
                    for i=1:num_oldBeadfiles
                        file=norm_cytof(fullfile(oldBead_dir,oldBead_files{i}),sbtrct);
                        file=file.find_bead_channels(bead_metals);
                        oldBead_vals{i}=file.data(:,file.bead_channels);
                        clear file
                    end
                    all_oldBeads=cat(1,oldBead_vals{:});
                else
                    file=norm_cytof(fullfile(oldBead_dir,oldBead_files),sbtrct);
                    file=file.find_bead_channels(bead_metals);
                    all_oldBeads=file.data(:,file.bead_channels);
                end
                bead_means=median(all_oldBeads);
            end
        else
            bead_means=median(allbeads(:,2:(num_beads+1)));
        end
        
    end

    function [normedbeads,smoothed,normedsmoothed]=normalizeFiles
        
        normedbeads=cell(1,num_files);
        smoothed=cell(1,num_files);
        normedsmoothed=cell(1,num_files);
        defAns={''};
        
        for i=1:num_files
            file=norm_cytof(filenames{i},sbtrct);
            file=file.find_bead_channels(bead_metals);
            file=file.assign_bead_data(bead_data{i}(:,end));
            file=file.smooth_beads(bead_window);
            file=file.compute_bead_slopes(bead_means);
            file=file.make_times_unique;
            file=file.correct_channels;  %changes values of file.data
            
            smoothed{i}=file.smoothed_beads; %these were smoothed before normalization
            file=file.assign_bead_data(bead_data{i}(:,end));  %moving normalized beads into bead_data property
            normedbeads{i}=file.bead_data;
            
            file=file.smooth_beads(bead_window);    %smoothing normalized beads
            normedsmoothed{i}=file.smoothed_beads;
            
            if rmv %&& sum(file.beads)>num_beads
                file=file.remove_beads(questpos,defAns);  %removes events within selected distance of beads
                
                if ~isempty(file.removed_events)
                    defAns={num2str(file.mahal_cutoff)};
                    fca_writefcs([bead_dir f{i}(1:end-4) '_normalized_removedEvents.fcs'],file.removed_events,file.markernames,file.channelnames,file.header);
                end
                fprintf(fid,'Removing %g events from %s with beadDist <= %g\n',size(file.removed_events,1),f{i},file.mahal_cutoff);
            end
            
            fca_writefcs([bead_dir f{i}(1:end-4) '_beads.fcs'],bead_data{i}(:,[1 end 2:end-1]),['Time' 'EventNum' file.channelnames(file.bead_channels)],['Time' 'EventNum' file.channelnames(file.bead_channels)],file.header);
            fca_writefcs([norm_dir f{i}(1:end-4) '_normalized.fcs'],file.data,file.markernames,file.channelnames,file.header);
            
            clear file
        end
        
    end

    function plotBeforeAndAfter
        %plot medians before and after (new version)
        
        finalfig=figure('color',[1 1 1],'NumberTitle','off','Name','Close Figure to Exit Program');
        ticks=load('axticks.mat');
        axticks=ticks.axticks;
        
        ax=zeros(1,2);
        ttl=zeros(1,2);
        leg=zeros(1,2);
        ax(1)=subplot(3,1,1);
        hold on
        ind=0;
        s=zeros(num_files,1);
        smoothedmeds=zeros(num_files,num_beads);
        normedmeds=zeros(num_files,num_beads);
        
        for i=1:num_files
            s(i)=size(smoothed{i},1);
            %     smoothedmeds(i,:)=median(smoothed{i});
            %     normedmeds(i,:)=median(normedsmoothed{i});
            smoothedmeds(i,:)=median(bead_data{i}(:,2:num_beads+1));
            normedmeds(i,:)=median(normedbeads{i}(:,2:num_beads+1));
            
        end
        if num_files>1
            plot(asinh(1/5*smoothedmeds),'o-','markersize',4,'linewidth',1)
            xl=[0.5 num_files+0.5];
        else
            s=size(smoothed{1},1);
            plot(1:s,asinh(1/5*smoothed{1}),'linewidth',1.5)
            if s>1
                xl=[1 s];
            else
                xl=[0 2];
            end
            % plot(ones(1,num_beads),asinh(1/5*smoothedmeds'),'o','markersize',4,'linewidth',1)
        end
        
        
        yl=get(gca,'ylim');
        
        filemeds=1:num_files;
        
        
        set(gca,'xlim',xl,'ylim',yl,'ytick',axticks.loc,'yticklabel',axticks.labels,'ygrid','on','xtick',[],'box','on','fontsize',10);
        % xlabel('time','fontsize',12)
        ttl(1)=title('Before Normalization','fontsize',12);
        leg(1)=legend(bead_metals,'fontsize',10,'location','eastoutside');
        
        ax(2)=subplot(3,1,2);
        hold on
        txt=zeros(1,num_files);
        if num_files>1
            plot(asinh(1/5*normedmeds),'o-','markersize',4,'linewidth',1)
        else
            s=size(normedsmoothed{1},1);
            plot(1:s,asinh(1/5*normedsmoothed{1}),'linewidth',1.5)
            %  plot(ones(1,num_beads),asinh(1/5*normedmeds'),'o','markersize',4,'linewidth',1)
        end
        for i=1:num_files
            p=strfind(f{i},'_cells_found');
            if isempty(p)
                p=length(f{i})+1;
            end
            txt(i)=text(filemeds(i),yl(1)-0.1*diff(yl),f{i}(1:p(1)-1),'rotation',45,'horizontalalignment','right','interpreter','none','fontsize',10);
        end
        set(gca,'xlim',xl,'ylim',yl,'ytick',axticks.loc,'yticklabel',axticks.labels,'ygrid','on','xtick',[],'box','on','fontsize',10);
        % xlabel('time','fontsize',12)
        ttl(2)=title('After Normalization','fontsize',12);
        leg(2)=legend(bead_metals,'fontsize',10,'location','eastoutside');
        
        set(finalfig,'paperposition',[0 0 6 3 ],'papersize',[6 3])
        
        set(ax,'fontsize',6)
        set(txt,'fontsize',6)
        set(leg,'fontsize',6)
        set(ttl,'fontsize',6)
        
        saveas(finalfig,[norm_dir filesep 'beads_before_and_after.png'],'png')
        
        set(ax,'fontsize',10)
        set(txt,'fontsize',10)
        set(leg,'fontsize',10)
        set(ttl,'fontsize',12)
        
        toofew=s<=1;
        smoothed(toofew)=[];
        normedsmoothed(toofew)=[];
        
        cat_before=cat(1,smoothed{:});
        max_before=max(cat_before);
        min_before=min(cat_before);
        
        cat_after=cat(1,normedsmoothed{:});
        max_after=max(cat_after);
        min_after=min(cat_after);
        
        r_before=mean(max_before./min_before);
        r_after=mean(max_after./min_after);
        
        fprintf('\n Bead Fractional Range Before = %0.5f\n',r_before);
        fprintf('\n Bead Fractional Range After = %0.5f\n',r_after);
        
        fprintf(fid,'Bead Fractional Range Before = %0.5f\n',r_before);
        fprintf(fid,'Bead Fractional Range After = %0.5f\n',r_after);
        fclose(fid);
       
    end

end


function [figpos,questpos] = setFigurePositions()
% set the default figure window to be full width of the screen and not too tall

nm=get(0,'monitorpositions');
if size(nm,1)==1 %one monitor, so screensize coordinates start at actual monitor corner
    ss=get(0,'screensize');
    figpos=[ss(1)+5 ss(4)/2 ss(3)-10 ss(4)/2.5];
    questpos=[ss(3)/2-200 ss(4)/4 400 200];
else
    defpos=get(0,'defaultfigureposition'); %screensize coords are inconsistent so give up on that
    figpos=[defpos(1)-defpos(3) defpos(2) 3*defpos(3) defpos(4)];
    questpos=[defpos(1) defpos(2)-defpos(4)/2 400 200];
end

end

function [data_dir, f,num_files]=getFileList
% load in list of filenames

start_dir='~';

data_dir = uigetdir(start_dir,'Choose Folder');
if data_dir==0
    return
else
    data_dir=[data_dir filesep];
    f=dir([data_dir '*.fcs']);
    f={f.name};
    num_files=length(f);
    if num_files == 0
        h=warndlg('No fcs files found!' ,'Empty folder!');
        waitfor(h)
        drawnow
        return
    end
end

end

function timeOrder=orderFilenames(filenames)
% put the filenames in time order (this will not be true time order if the files have been re-extracted)

num_files=length(filenames);
h=cell(1,num_files);
for i=1:num_files
    h{i}=get_fcs_hdr(filenames{i});
end
timeOrder=get_fcs_order(h);

end

function cellData=concatenateTimestamps(cellData)

for i=2:length(cellData)
    cellData{i}(:,1)=cellData{i}(:,1)-cellData{i}(1,1)+mean(diff(cellData{i}(:,1)))+cellData{i-1}(end,1);
end

end

