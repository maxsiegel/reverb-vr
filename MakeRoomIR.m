
clear all; close all; clc

%% Dash
Rm=[40 40 8]; %Room dimensions [LR FB UP]
sx=[30 30 2]; % Source position
lx=[5 5 2]; % Listener position

wflt=-10; % Reflection co-efficient (this should be a filter -peaked between
          % 200-2000Hz)

fs=44100;
c=340;

%% Crunch

% scroll through number of reflections
h=zeros(2*fs,1);
for jrf=0:6;
    %% move the source accordingly
    % compute the difference vectors by which the image source can move
    Df1=2*lx;
    Df2=2*(Rm-lx);
    % pre-allocate
    prms=zeros(1,3);
    % scroll through the number of ways the number of reflections can be
    % split among 3 wall pairs
    for jrf2=1:jrf;
        tprms=zeros(1,3);
        tprms(1)=jrf;
        el=1;
        while(tprms(1))>(jrf-jrf2+1); el=el+1; el=rem(el,3)+1;
          tprms(1)=tprms(1)-1;
          tprms(el)=tprms(el)+1;
        end
        tprms=perms(tprms);
        % remove repetitions
        cnt=0;
        while cnt<size(tprms,1)-1; cnt=cnt+1;
            cnt2=cnt;
            while cnt2<size(tprms,1); cnt2=cnt2+1;
                if tprms(cnt,:)==tprms(cnt2,:);
                    tprms(cnt2,:)=[];
                    cnt2=cnt2-1;
                end
            end
        end
        prms=[prms; tprms];
    end
    if jrf~=0; prms(1,:)=[]; end
     % scroll through the permutations of wall reflections
     for jprm=1:size(prms);
         % there are two sides for every-pairing
         if jrf==0; Nrf=1; else Nrf=2; end
         for jdr=1:Nrf
            tprms=prms(jprm,:);
            if jdr==1;
                tDf=-rem(tprms,2).*Df1+abs(rem(tprms,2)-1).*Df2;
            elseif jdr==2;
                tprms=-tprms;
                tDf=rem(tprms,2).*Df2-abs(rem(tprms,2)-1).*Df1;
            end
            tDf=tDf.*tprms;
            tprms2=max(tprms-1,0);
            tlx=lx+tprms2.*Rm+tDf;
            % find path length
            L=rms(sx-tlx);
            % find the direct arrival
            if jrf==0; L0=L; end
            if (jrf==0||L>L0)
              % find the arrival time
              t=L/c*fs;
              % add to the time series
              h(ceil(t))=h(ceil(t))+10^(jrf*wflt/20)/(L^2);
              % plot
              figure(101); plot(h); axis tight; set(gca,'xlim',[-10 25e3]);
              title(sprintf('%d reflections',jrf)); drawnow;
              fprintf(['No of reflections=%d, walls=[%d %d %d], dst=[%2.1f ' ...
                       '%2.1f %2.1f], tlx=[%2.1f %2.1f %2.1f], rms=%2.1f\n'],jrf,tprms(1),tprms(2),tprms(3),tDf(1),tDf(2),tDf(3),tlx(1),tlx(2),tlx(3),L);
            end
         end
     end
end
h=h/max(abs(h));

% convolve with itself
ndx=find(h>0);
th=h(min(ndx):max(ndx));
nh=h;
for jrp=1:20
    th=th(2:end);
    th=th(randperm(length(th)));
    th(find(th~=0))=th(find(th~=0)).*randn(size(find(th~=0)));
    th=[1; th];
    nh=RIRcnv(nh,fs,th,fs,1);
    figure(102); subplot(2,1,1); plot(nh); axis tight; 
    subplot(2,1,2); plot(20*log10(abs(nh))); axis tight; 
    set(gca,'ylim',[-100 0])
    title(sprintf('%d convolutions',jrp)); drawnow; 
end
