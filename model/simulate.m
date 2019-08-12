% run simulation and print output
clear; close all;

addpath('/Users/danielborrus/Documents/GitHub/wavelets/wave_matlab/')

%% Parameters

P.tmax = 12;   % maximum simulation time in secconds
P.dt   = 0.001; % time step

P.N    = 10;   % number of neurons
P.p    = 1;    % period of spiking
P.vp   = 2;    % STD of spiking
P.buf  = 2;    % amount of buffer time before/after spiking

%% Presim calculations and initializations

P.nsteps  = P.tmax/P.dt;
P.nspikes = P.tmax/P.p - (P.buf/P.p-1) - (P.buf/P.p);

if  2*P.vp>=P.buf
    warning(['A period variance (P.vp) larger than or close too the front buffer ' ...
             '(P.buf) may push a spike before the simulation start, causing an error. '...
             ' Recommend P.vp<<P.buf'])
end

%% Blocks

[sumt,tra] = GenSpikTrain(P);

[ctrain,fil,Out] = Convolut_custom(sumt,P);

[cfs,frq] = Wanal(ctrain,P);

%% Figure

if 1
    
figure('Position',[100 700 1000 700])
subplot(2,1,1); hold on;
plot(P.dt:P.dt:length(ctrain)*P.dt,ctrain,'linewidth',2);
title('Synaptic Drive')
ylabel('')
axis tight
ax1 = gca;
ax1.Position=ax1.Position - [0 0 .1 0];
dim = [0.82 0.6 0.3 0.3];
str = {[mat2str(P.N) ' neurons'],['Period = ' mat2str(P.p) ' sec'],['St. Dev. = ' mat2str(P.vp) ' sec']};
annotation('textbox',dim,'String',str,'FitBoxToText','on','edgecolor','k');

subplot(2,1,2); hold on;
Fs = 1/P.dt;
tms = (0:numel(ctrain)-1)/Fs;
surface(tms,frq,abs(cfs))
title('Wavelet Scalogram')
axis tight
shading flat
xlabel('Time (s)')
ylabel('Frequency (Hz)')
ylim([4 64])
set(gca,'yscale','log')
yticks([4 16 64])
ax2 = gca;
ax2.Position=ax2.Position - [0 0 .1 0];

Ytick = strsplit(num2str(get(gca, 'Ytick')));
set(gca,'YtickLabel', Ytick);

cb=colorbar('location','EastOutside');
cb.Position=cb.Position + [.1 0 0 0];

set(findall(gcf,'-property','FontSize'),'FontSize',20)

saveas(gcf,'wavalysis','png'); close


figure('Position',[100 700 1000 700])
subplot(4,1,1); hold on;
plot(tra); title('Individual spike trains')
yticks(1)
yticklabels('spike'); axis tight
subplot(4,1,2); hold on;
plot(sumt); title('Summated spike train'); axis tight
subplot(4,1,3); hold on;
plot(fil); title(['EPSP shape (n=' mat2str(Out.n) ', tau=' mat2str(Out.tau) ')']); axis tight
subplot(4,1,4); hold on;
plot(ctrain); title('Convoluted spike train'); axis tight
set(findall(gcf,'-property','FontSize'),'FontSize',16)

saveas(gcf,'extra','png'); close;

figure
cwt(ctrain,'amor',1/P.dt);
set(findall(gcf,'-property','FontSize'),'FontSize',16)
saveas(gcf,'rawwavelet','png'); close;


end


%% Functions

function [sumtrains,trains] = GenSpikTrain(P)
% generate a spike train trace

% parameters
dofig = 0;

% precalc stuff
P.tmax_raw = P.tmax/P.dt;
P.dt_raw   = P.dt/P.dt;
P.p_raw    = P.p/P.dt;
P.vp_raw   = P.vp/P.dt;
P.buf_raw  = P.buf/P.dt;

% main chunk

% set unjittered (raw) spike times
spike_ts = (P.p_raw * ((P.buf/P.p):P.tmax/P.p-(P.buf/P.p)) )';
% increase number of columns from 1 --> N
spike_ts = repmat(spike_ts, [1,P.N] );
% I realized later this line is unneccesary thanks to the matrix addition
% next, but i keep it in for easier to read code


% add jitter to (raw) spike time
spike_ts = spike_ts + P.vp_raw * randn(P.nspikes,P.N);
% round to integers for indexing
spike_ts = round(spike_ts);


% empty array of neurons spike train (0 - no spike, 1 - spike)
trains = zeros(P.nsteps,P.N);

% fill the array
for i = 1:P.N; trains(spike_ts(:,i),i) = 1; end

% output
sumtrains = sum(trains,2);

if dofig
   figure
   subplot(2,1,1); hold on; title('Single neuron output')
   plot(P.dt:P.dt:length(sumtrains)*P.dt,trains(:,1))
   subplot(2,1,2); hold on; title('Summed output of all neurons')
   plot(P.dt:P.dt:length(sumtrains)*P.dt,sumtrains)
end

end

function [ctrain,fil,Out] = Convolut_custom(sumtrains,P)
% convolute the summed spikes train with a nice decay model


%fil = (exp(linspace(0.6931,0,P.dcay/P.dt/1000))-1);
fil_thefunction = @(n,tau,t) (1/(factorial(n)*tau))*((t/tau).^n).*exp(-t/tau);

n = 1;tau = 20;
fil = fil_thefunction(n,tau,0:180);

ctrain = conv(sumtrains,fil);


if 0
figure
subplot(2,1,1)
plot(fil)
subplot(2,1,2)
plot(ctrain)
end

Out.tau = tau;
Out.n = n;

end

function [cfs,frq] = Wanal(ctrain,P)

% figure('Position',[100 700 1000 700])
% subplot(2,1,1); hold on;
% plot(P.dt:P.dt:length(ctrain)*P.dt,ctrain);
% title('Network Activity')
% ylabel('')
% axis tight
% ax1 = gca;
% ax1.Position=ax1.Position - [0 0 .1 0];
% dim = [0.82 0.6 0.3 0.3];
% str = {[mat2str(P.N) ' neurons'],['Period = ' mat2str(P.p) ' sec'],['St. Dev. = ' mat2str(P.vp) ' sec'],['Decay = ' mat2str(P.dcay) ' msec']};
% annotation('textbox',dim,'String',str,'FitBoxToText','on','edgecolor','k');

% subplot(2,1,2); hold on;
[cfs,frq]= cwt(ctrain,'amor',1/P.dt);
%cwt(ctrain,'morse',1000);
% 
% Fs = 1000;
% tms = (0:numel(ctrain)-1)/Fs;
% %subplot(2,1,2)
% surface(tms,frq,abs(cfs))
% title('Wavelet Scalogram')
% axis tight
% shading flat
% xlabel('Time (s)')
% ylabel('Frequency (Hz)')
% ylim([4 64])
% set(gca,'yscale','log')
% yticks([4 16 64])
% ax2 = gca;
% ax2.Position=ax2.Position - [0 0 .1 0];
% 
% Ytick = strsplit(num2str(get(gca, 'Ytick')));
% set(gca,'YtickLabel', Ytick);


% 

% 
% cb=colorbar('location','EastOutside');
% cb.Position=cb.Position + [.1 0 0 0];

% 
% set(findall(gcf,'-property','FontSize'),'FontSize',20)
%wavelet(Y,DT,PAD,DJ,S0,J1,MOTHER,PARAM)
%WAVE = wavelet(ctrain,.001);



end
