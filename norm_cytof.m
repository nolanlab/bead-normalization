classdef norm_cytof
    
    properties
        data
        header
        channelnames
        markernames
        time_col
        cell_length_col
        dna_channels
        marker_channels
        
        bead_metals
        bead_channels
        beads
        xp
        yp
        
        bead_data
        smoothed_beads
        smoothed_time
        bead_slopes
        
        ticks
        ticklabels
        
        numcm
        removed_events
        mahal_cutoff
    end
    
    methods
        
        function obj = norm_cytof(filename,sbtrct)
            
            
            
            [obj.data, obj.header] = fca_readfcs(filename);
            if sbtrct
                obj.data(:,3:end)=obj.data(:,3:end)-100;
            end
            
            obj.channelnames={obj.header.par.name};
            obj.markernames={obj.header.par.name2};
            
            % find time, cell length and DNA channels
            obj.time_col = find(strcmpi('time',obj.channelnames));
            obj.cell_length_col = find(strcmpi('cell_length',obj.channelnames));
            
            dna_1=find(~cellfun(@isempty,strfind(obj.channelnames,'Ir191')));
            dna_2=find(~cellfun(@isempty,strfind(obj.channelnames,'Ir193')));
            
            obj.dna_channels=[dna_1 dna_2];
            
            obj.marker_channels=1:length(obj.channelnames);
            obj.marker_channels([obj.time_col obj.cell_length_col])=[];
            
            %load in axes ticks and labels
            load('axticks.mat')
            obj.ticks=axticks.loc;
            obj.ticklabels=axticks.labels;
        end
        
        function obj = find_bead_masses(obj)
            
            global num_beads
            
            m=obj.channelnames;
            
            default_masses={'La139','Pr141','Tb159','Tm169','Lu175'};
            num_defaults=length(default_masses);
            bead_cols=zeros(1,num_defaults);
            
            for j=1:num_defaults
                bci=find(~cellfun(@isempty,strfind(m,default_masses{j})));
                if ~isempty(bci)
                    bead_cols(j)=bci;
                end
            end
            
            massfig=figure('KeyPressFcn',@mass_return);
            col_list=uicontrol('parent',massfig,...
                'style','listbox',...
                'string',m,...
                'units','normalized',...
                'position',[0.1 0.3 0.8 0.6],...
                'max',2,...
                'min',0,...
                'value',bead_cols(bead_cols ~= 0)); %it's possible that not all default bead masses were found
            uicontrol('parent',massfig',...
                'style','pushbutton',...
                'string','Choose Bead Masses',...
                'units','normalized',...
                'position',[0.4 0.1 0.2 0.1],...
                'Callback',{@set_bead_metals})
            
            waitfor(massfig)
            
            function mass_return(src,event)
                
                if strcmp(event.Key,'return')
                    set_bead_metals
                end
            end
            
            function set_bead_metals(varargin)
                obj.bead_channels=get(col_list,'value');
                bead_colnames=obj.channelnames(obj.bead_channels);
                num_beads=length(bead_colnames);
                obj.bead_metals=cell(1,num_beads);
                for i=1:num_beads
                    metal=obj.channelnames{obj.bead_channels(i)};
%                     p1=strfind(metal,'(');
%                     p2=strfind(metal,')');
%                     obj.bead_metals{i}=metal(p1(end):p2(end));
                    obj.bead_metals{i}=metal;
                end
                close(massfig);
            end
            
        end
        
        function obj = find_bead_channels(obj,bead_metals)
            
            m=obj.channelnames;
            num_beads=length(bead_metals);
            bead_cols=zeros(1,num_beads);
            for i=1:num_beads
                bci=find(~cellfun(@isempty,strfind(m,bead_metals{i})));
                if ~isempty(bci)
                    bead_cols(i)=bci;
                end
            end
            obj.bead_channels=bead_cols;
        end
        
        function obj = manual_bead_identify_old(obj,varargin)
            
            dnadata=asinh(1/5*(obj.data(:,obj.dna_channels(1))));
            beaddata=asinh(1/5*(obj.data(:,obj.bead_channels)));
            
            beads=true(size(dnadata));
%             fig=figure('position',[27 290 1230 341]);

            ax=zeros(1,5);
            for i=1:5
                ax(i);subplot(1,5,i);
                plot(ax(i),beaddata(:,i),dnadata,'.','markersize',2)
            end
            
            drawnow;
            
            if nargin>1
                gates=norm_cytof.gate_beads(fig,ax,varargin{1},varargin{2});
            else
                gates=norm_cytof.gate_beads(fig,ax);
            end
            
            pause
            
            xp=cell(1,5);
            yp=cell(1,5);
            for i=1:5
                xp{i}=get(gates(i),'xdata');
                yp{i}=get(gates(i),'ydata');
                
                beads=beads & inpolygon(beaddata(:,i),dnadata,xp{i}',yp{i}');
            end
            
            obj.xp=xp;
            obj.yp=yp;
            obj.beads= beads;
            close(fig)
            drawnow;
        end
        
        function obj = manual_bead_identify(obj,xp,yp)
            
          
            
            global figpos
            global num_beads
            
            
            dnadata=asinh(1/5*(obj.data(:,obj.dna_channels(1))));
            beaddata=asinh(1/5*(obj.data(:,obj.bead_channels)));
            
            %if too many cells, sample for visualization
            num_cells=size(obj.data,1);
            if num_cells > 100000
                inds=randsample(num_cells,100000);
                sampled_dnadata=dnadata(inds,:);
                sampled_beaddata=beaddata(inds,:);
            else      
                sampled_dnadata=dnadata;
                sampled_beaddata=beaddata;
            end
            
            
            fig=figure('position',figpos,...
                'KeyPressFcn',@replot,...
                'ResizeFcn',@change_figpos,...
                'NumberTitle','off',...
                'Name',['Press Enter to Confirm Gates (Optionally, Press Space Bar First to View Selected Cells): ' obj.header.filename]);
            ax=zeros(1,num_beads);
            ph=zeros(1,num_beads);
            for i=1:num_beads
                ax(i)=subplot(1,num_beads,i);
                hold on
                ph(i)=plot(sampled_beaddata(:,i),sampled_dnadata,'.','color',[0.3 0.3 .3],'markersize',2);
                set(ax(i),'xtick',obj.ticks,'xticklabel',obj.ticklabels,'ytick',obj.ticks,'yticklabel',obj.ticklabels)
                set(ax(i),'xtick',obj.ticks,'xticklabel',[],'ytick',obj.ticks,'yticklabel',[])
                xlabel(obj.channelnames{obj.bead_channels(i)})
                ylabel(obj.channelnames{obj.dna_channels(1)});
            end
            
            
            if exist('xp','var') && ~isempty(xp)
                gates=norm_cytof.gate_beads(fig,ax,xp,yp);
            else
                gates=norm_cytof.gate_beads(fig,ax);
            end
            
            waitfor(fig)
            
            
            function replot(src,event)
                if exist('gates','var')  %these won't exist yet if you hit a key too soon
                
                
                if strcmp(event.Key,'space')
                    beadcolor=[1 0 0];
                    
                    beads=true(size(sampled_dnadata));
                    xp=cell(1,num_beads);
                    yp=cell(1,num_beads);
                    for j=1:num_beads
                        xp{j}=get(gates(j),'xdata');
                        yp{j}=get(gates(j),'ydata');
                        
                        beads = beads & inpolygon(sampled_beaddata(:,j),sampled_dnadata,xp{j},yp{j});
                    end
                    
                    for j=1:num_beads
                        ch=findall(ax(j),'color',beadcolor);
                        if ~isempty(ch)
                            set(ch,'xdata',sampled_beaddata(beads,j),'ydata',sampled_dnadata(beads))
                            set(ax(j),'children',[gates(j) ch ph(j)]);
                        else
                            ch=plot(ax(j),sampled_beaddata(beads,j),sampled_dnadata(beads),'+','color',beadcolor,'markersize',4);
                            set(ax(j),'children',[gates(j) ch ph(j)]);
                        end
                    end
                    drawnow;
                elseif strcmp(event.Key,'return')
                    
                    %use all, rather than sampled, data now
                    beads=true(size(dnadata));
                    xp=cell(1,num_beads);
                    yp=cell(1,num_beads);
                    for j=1:num_beads
                        xp{j}=get(gates(j),'xdata');
                        yp{j}=get(gates(j),'ydata');
                        
                        beads = beads & inpolygon(beaddata(:,j),dnadata,xp{j},yp{j});
                    end
                    
                    obj.xp=xp;
                    obj.yp=yp;
                    obj.beads= beads;
                    close(fig)
                    drawnow;
                end
                else
                    return
                end
            end
            
            function change_figpos(src,event)
                figpos=get(fig,'position');
            end
            
        end
        
        function  obj = add_bead_col(obj)
            
            y=zeros(size(obj.data,1),size(obj.data,2)+1);
            y(:,1:end-1)=obj.data;
            y(:,end)=obj.beads;
            obj.data=y;
            obj.channelnames=[obj.channelnames 'beads'];
            obj.markernames=[obj.markernames 'beads'];
        end
        
        function obj = assign_bead_data(obj,bead_eventnums)
            
            if nargin>1
                obj.bead_data=obj.data(bead_eventnums,[1 obj.bead_channels]);
                obj.beads=false(size(obj.data,1),1);
                obj.beads(bead_eventnums)=true;
            elseif any(obj.beads)
                obj.bead_data=obj.data(obj.beads,[1 obj.bead_channels]);
            else
                error('bead eventnums or obj.beads not assigned')
            end
        end
        
        function obj = smooth_beads(obj, window_size)
            
            global num_beads
            
            num_events = size(obj.bead_data, 1);
            
            if num_events >= window_size
                smoothed_beads = zeros(num_events-window_size, num_beads);
                smoothed_time = zeros(num_events-window_size, 1);
                
                for i=1:num_events-window_size
                    
                    smoothed_beads(i,:) = median(obj.bead_data(i:i+window_size, 2:(num_beads+1)));
                    smoothed_time(i) = median(obj.bead_data(i:i+window_size, 1));
                    
                end
                
            else
                
                smoothed_beads=median(obj.bead_data(:,2:(num_beads+1)));
                smoothed_time = median(obj.bead_data(:,1));
                
            end
            
            obj.smoothed_beads=smoothed_beads;
            obj.smoothed_time=smoothed_time;
            
            
        end
        
        function obj = compute_bead_slopes(obj,baseline)
            
            global num_beads
            
            if nargin<2
                baseline=mean(obj.bead_data(2:(num_beads+1)));
            end
            
            s=size(obj.smoothed_beads,1);
            y=repmat(baseline,[s 1]);
            
            obj.bead_slopes=sum(obj.smoothed_beads.*y,2)./sum(obj.smoothed_beads.^2,2);
            
        end
        
        
        function obj = make_times_unique(obj)
            
            t=obj.smoothed_time;
            dt=diff(t);
            repeats=find(dt==0);
            
            if ~isempty(repeats)
                ri=0;
                
                if repeats(end)< (length(t)-1)
                    for i=1:length(repeats)
                        if repeats(i)==ri
                            continue
                        end
                        ri=repeats(i);
                        next_t=find(t>t(ri),1,'first');
                        new_ts=linspace(t(ri),t(next_t),next_t-ri+1);
                        t(ri:next_t)=new_ts;
                    end
                else  %special case because can't go beyond end because next_t would be empty
                    L=find(diff(t(repeats)),1,'last');
                    for i=1:L
                        if repeats(i)==ri
                            continue
                        end
                        ri=repeats(i);
                        next_t=find(t>t(ri),1,'first');
                        new_ts=linspace(t(ri),t(next_t),next_t-ri+1);
                        t(ri:next_t)=new_ts;
                    end
                    for i=repeats(L+1)+1:length(t)
                        t(i)=t(i-1)+t(i-1)-t(i-2);
                    end
                end
                
                obj.smoothed_time=t;
                
            end
        end
        
        
        function obj=correct_channels(obj)
            
            if length(obj.smoothed_time)>1
                
                ir=interp1(obj.smoothed_time, obj.bead_slopes, obj.data(:,obj.time_col),'linear');
                
                tstart2=find(~isnan(ir),1,'first');
                ir(1:tstart2)=ir(tstart2);
                
                tend=find(~isnan(ir),1,'last');
                ir(tend:end)=ir(tend);
                
            else
                
                ir=obj.bead_slopes*ones(size(obj.data(:,1)));
                
            end
            
            for j=obj.marker_channels  %CHANGING DATA VALUES FOR GOOD HERE
                obj.data(:,j)=obj.data(:,j).*ir;
            end
            
        end
        
        
        function obj = remove_beads(obj,questpos,defAns)
            
            global figpos
            global num_beads
            
            tdata=asinh(1/5*obj.data(:,obj.bead_channels));          
            
            if sum(obj.beads)>num_beads
                mdists=mahal(tdata,tdata(obj.beads,:));
                mdists=sqrt(mdists);
                
                %sample for visualization if too many cells
                num_cells=size(tdata,1);
                if num_cells>100000
                    inds=randsample(num_cells,100000);
                    tdata_sampled=tdata(inds,:);
                    mdists_sampled=mdists(inds,:);
                else
                    tdata_sampled=tdata;
                    mdists_sampled=mdists;
                end
            else
                
                %                obj.removed_events=[obj.data(~notbeads,:) mdists(~notbeads)];
                obj.data=[obj.data zeros(size(obj.data,1),1)];
                obj.channelnames=[obj.channelnames 'beadDist'];
                obj.markernames=[obj.markernames 'beadDist'];
                obj.mahal_cutoff=0;
                h=warndlg(['Too few beads to estimate distances in ' obj.header.filename ': you must remove them manually.'] ,'Sparse beads!');
                waitfor(h)
                drawnow
                return
            end
            
            
            sc=get(0,'screensize');
            fig=figure('position',figpos,...
                'KeyPressFcn',@mahal_cutoff,...
                'ResizeFcn',@change_figpos,...
                'NumberTitle','off',...
                'Name',obj.header.filename);
            %             fig=figure('position',[1 0.3*sc(4) sc(3) sc(3)/3]);
            
            nbins=20;
            [~,bin]=histc(mdists_sampled,0:nbins);
            
            cm=flipud(jet(nbins+1));
            colormap(cm);
            obj.numcm=nbins+1;
            
            
            if num_beads>2
                
                combos=combntns(1:num_beads,2);
                
                ax=zeros(1,3);
                ax(1)=axes('position',[0.0625 0.1 0.25 0.75]);
                hold on
                plot(tdata_sampled(bin==0,combos(1,1)),tdata_sampled(bin==0,combos(1,2)),'+','color',cm(end,:),'markersize',2)
                for j=1:nbins
                    i=nbins+1-j;
                    plot(tdata_sampled(bin==i,combos(1,1)),tdata_sampled(bin==i,combos(1,2)),'+','color',cm(i,:),'markersize',2)
                end
                
                xlabel(obj.channelnames{obj.bead_channels(combos(1,1))})
                ylabel(obj.channelnames{obj.bead_channels(combos(1,2))})
                
                ax(2)=axes('position',[0.375 0.1 0.25 0.75]);
                hold on
                plot(tdata_sampled(bin==0,combos(end-1,1)),tdata_sampled(bin==0,combos(end-1,2)),'+','color',cm(end,:),'markersize',2)
                for j=1:nbins
                    i=nbins+1-j;
                    plot(tdata_sampled(bin==i,combos(end-1,1)),tdata_sampled(bin==i,combos(end-1,2)),'+','color',cm(i,:),'markersize',2)
                end
                
                xlabel(obj.channelnames{obj.bead_channels(combos(end-1,1))})
                ylabel(obj.channelnames{obj.bead_channels(combos(end-1,2))})
                
                ax(3)=axes('position',[0.6875 0.1 0.25 0.75]);
                hold on
                plot(tdata_sampled(bin==0,combos(end,1)),tdata_sampled(bin==0,combos(end,2)),'+','color',cm(end,:),'markersize',2)
                for j=1:nbins
                    i=nbins+1-j;
                    plot(tdata_sampled(bin==i,combos(end,1)),tdata_sampled(bin==i,combos(end,2)),'+','color',cm(i,:),'markersize',2)
                end
                
                xlabel(obj.channelnames{obj.bead_channels(combos(end,1))})
                ylabel(obj.channelnames{obj.bead_channels(combos(end,2))})
                
                for i=1:3
                    set(ax(i),'xtick',obj.ticks,'xticklabel',obj.ticklabels,'ytick',obj.ticks,'yticklabel',obj.ticklabels)
                end
                
            else
                
                
                ax=axes('position',[0.375 0.1 0.25 0.75]);
                hold on
                for i=1:nbins
                    plot(tdata_sampled(bin==i,1),tdata_sampled(bin==i,2),'+','color',cm(i,:),'markersize',2)
                end
                plot(tdata_sampled(bin==0,1),tdata_sampled(bin==0,2),'+','color',cm(end,:),'markersize',2)
                xlabel(obj.channelnames{obj.bead_channels(1)})
                ylabel(obj.channelnames{obj.bead_channels(2)})
                
                
                set(ax,'xtick',obj.ticks,'xticklabel',obj.ticklabels,'ytick',obj.ticks,'yticklabel',obj.ticklabels)
            end
            
            
            cbax=axes('position',[0 0 1 1],'clim',[0 nbins+1],'visible','off');
            cb=colorbar('peer',cbax,'location','north','position',[0.0625 0.9 0.875 0.05]);
            set(get(cb,'title'),'string','Distance from Identified Beads')
            
            txt=['Type in the cutoff, then press Enter.\n' ...
                '(You can enter 0 if you do not want to apply a cutoff.)\n'];
            
            drawnow;
            %             if ~isdeployed
            %             commandwindow;
            %
            %             while (1)
            %                 try cutoff=input(txt);
            %                 catch exception
            %                     if ~isempty(exception.identifier)
            %                         sprintf('%s','Try again: Cutoff must be a number > 0.')
            %                     end
            %                 end
            %
            %                 if 0<= cutoff < Inf
            %                     notbeads=mdists>cutoff;
            %                     break
            %                 end
            %             end
            %             else
            QuestOpt.Position=questpos;
            if nargin<3
                defAns={''};
            end
            while (1)
                cutoff=str2double(myInputdlg({'Type in the cutoff, then press Enter. (You can enter 0 if you do not want to apply a cutoff.)'},...
                    'Choose Cutoff',1,defAns,QuestOpt));
                if 0<=cutoff && cutoff<Inf && ~isnan(cutoff)
                    break
                end
            end
            
            
            close(fig)
            drawnow;
            
            notbeads=mdists>cutoff;
            %             end
            obj.removed_events=[obj.data(~notbeads,:) mdists(~notbeads)];
            obj.data=[obj.data(notbeads,:) mdists(notbeads)];
            obj.channelnames=[obj.channelnames 'beadDist'];
            obj.markernames=[obj.markernames 'beadDist'];
            
            obj.mahal_cutoff=cutoff;
            
            function mahal_cutoff(src,event)
                shift=false;
                shift_size=5;
                switch event.Key
                    case 'uparrow'
                        shift=true;
                        obj.numcm=obj.numcm+shift_size;
                        
                    case 'downarrow'
                        shift=true;
                        obj.numcm=obj.numcm-shift_size;
                        
                end
                
                
                if shift && obj.numcm>1
                    ncm=obj.numcm-1;
                    [~,bin]=histc(mdists,0:ncm);
                    cm=flipud(jet(obj.numcm));
                    
                    cla(ax(1))
                    for i=1:ncm
                        plot(ax(1),tdata(bin==i,combos(1,1)),tdata(bin==i,combos(1,2)),'+','color',cm(i,:),'markersize',2)
                    end
                    plot(ax(1),tdata(bin==0,combos(1,1)),tdata(bin==0,combos(1,2)),'+','color',cm(end,:),'markersize',2)
                    
                    cla(ax(2))
                    for i=1:ncm
                        plot(ax(2),tdata(bin==i,combos(end-1,1)),tdata(bin==i,combos(end-1,2)),'+','color',cm(i,:),'markersize',2)
                    end
                    plot(ax(2),tdata(bin==0,combos(end-1,1)),tdata(bin==0,combos(end-1,2)),'+','color',cm(end,:),'markersize',2)
                    
                    cla(ax(3))
                    for i=1:ncm
                        plot(ax(3),tdata(bin==i,combos(end,1)),tdata(bin==i,combos(end,2)),'+','color',cm(i,:),'markersize',2)
                    end
                    plot(ax(3),tdata(bin==0,combos(end,1)),tdata(bin==0,combos(end,2)),'+','color',cm(end,:),'markersize',2)
                    
                    set(cbax,'clim',[0 obj.numcm])
                    drawnow;
                end
            end
            
            
            function change_figpos(src,event)
                figpos=get(fig,'position');
            end
            
        end
        
        
    end
    
    
    methods(Static)
        
        
        function gate=gate_beads(fig,ax,xp,yp)
            
            global num_beads
            
            
            if nargin<3
                
                xp=cell(1,num_beads);
                yp=cell(1,num_beads);
                
                xp{1}=[2 4 4 2 2];
                xp{2}=[2 4 4 2 2];
                yp{1}=[-1 -1 0.5 0.5 -1];
                yp{2}=[-1 -1 0.5 0.5 -1];
                
                
                if num_beads>2
                    for i=3:num_beads
                        xp{i}=[4 6 6 4 4];
                        yp{i}=[-1 -1 0.5 0.5 -1];
                    end
                end
                
                
                
            end
    
            
            gate=zeros(1,num_beads);
            for i=1:num_beads
                set(gcf,'currentaxes',ax(i));
                gate(i)=line(xp{i},yp{i},'marker', '.','color',[0 0.25 1],'buttondownfcn',{@start_drag,i});
%                 set(ax(i),'ButtonDownFcn',{@rect_drag,i});
            end
            
            drawnow;
%             
%             function rect_drag(src,event,axi)
%                 xd=get(gate(axi),'xdata')
%                 yd=get(gate(axi),'ydata');
%                 cp=get(ax(axi),'currentpoint');
%                 if inpolygon(cp(1,1),cp(1,2),xd,yd)
%                     display('yup!')
% %                 else
% %                     return
%                 end
%             end
%             
            
            function start_drag(src,event,axi)
                
                
                cp=get(ax(axi),'currentpoint');
%                 [~,pti]=min((xp{axi}-cp(1,1)).^2 + (yp{axi}-cp(1,2)).^2);  %which corner is selected
                
                xd=get(gate(axi),'xdata');
                yd=get(gate(axi),'ydata');
                
                [~,pti]=min((xd-cp(1,1)).^2 + (yd-cp(1,2)).^2);  %which corner is selected

                set(fig,'WindowButtonMotionFcn',{@drag,axi,pti,xd,yd},'WindowButtonUpFcn',@stopdrag)
                
            end
            
            function drag(src,event,axi,pti,xd,yd)
                if pti==1
                    cp=get(ax(axi),'currentpoint');
                    
                    xd([1 4 5])=cp(1,1);
                    yd([1 2 5])=cp(1,2);
                    set(gate(axi),'xdata',xd,'ydata',yd);
                elseif pti==2
                    cp=get(ax(axi),'currentpoint');
                    yd=get(gate(axi),'ydata');
                    xd([2 3])=cp(1,1);
                    yd([1 2 5])=cp(1,2);
                    set(gate(axi),'xdata',xd,'ydata',yd);
                elseif pti==3
                    cp=get(ax(axi),'currentpoint');
                    
                    xd([2 3])=cp(1,1);
                    yd([3 4])=cp(1,2);
                    set(gate(axi),'xdata',xd,'ydata',yd);
                elseif pti==4
                    cp=get(ax(axi),'currentpoint');
                    
                    xd([1 4 5])=cp(1,1);
                    yd([3 4])=cp(1,2);
                    set(gate(axi),'xdata',xd,'ydata',yd);
                end
            end
            
            function stopdrag(varargin)
                set(fig,'WindowButtonMotionFcn',[])
            end
            
        end
    end
    
end













