% run simulation and print output
clear; close all;

%% Parameters

P.tmax = 20;   % maximum simulation time in secconds
P.dt   = 0.001; % time step
P.t    = 0:P.dt:P.tmax-P.dt; % time array
P.tct  = [0:P.dt:P.tmax-P.dt+(180*P.dt)]'; % time array for convoluted 
                                           % signal (it's longer due to 
                                           % the convolution).

P.N    = 10;   % number of neurons
P.p    = 1;    % period of spiking
P.vp   = 0.2;    % STD of spiking
P.buf  = 5;    % amount of buffer time before/after spiking

%% Presim calculations and initializations

P.nsteps  = P.tmax/P.dt;
P.nspikes = P.tmax/P.p - (P.buf/P.p-1) - (P.buf/P.p);

if  2*P.vp>=P.buf
    warning(['A period variance (P.vp) larger than or close too the front buffer ' ...
             '(P.buf) may push a spike before the simulation start, causing an error. '...
             ' Recommend P.vp<<P.buf'])
end

%% Blocks

% generate both spike trains
[sumt1,tra1] = GenSpikTrain(P);
[sumt2,tra2] = GenSpikTrain(P);

% convule the spike trains an alpha function
[ctrain1,fil,Out] = Convolut_custom(sumt1,P);
[ctrain2,~,~] = Convolut_custom(sumt2,P);

% Do the cross correlation
[corrCo,lags,trimdsignals] = CrossThem(ctrain1,ctrain2,P);

%% Figure

if 1
    
    figure('Position',[100 100 1000 700])
    subplot(3,3,[1 3])
    plot(P.t,sumt1+max(sumt2)+1,'b',P.t,sumt2,'r','linewidth',2);
    title('Two spike trains')
    legend('Binary input to neuron 1','Binary input to neuron 2')
    xlim([0 P.tmax])
    subplot(3,3,5)
    plot(fil,'k');
    title('Alpha function, or "The perfect EPSP"')
    set(gca,'XTick',[], 'YTick', []);
    subplot(3,3,[7 9])
    plot(P.tct,ctrain1+max(ctrain2)*1.2,'b',P.tct,ctrain2,'r','linewidth',2);
    title('Convoluted signal, i.e. the synaptic input')
    legend('Input to neuron 1','Input to neuron 2')
    xlim([0 P.tmax]) 
    saveas(gcf,'introduction','png'); %close
    
    
    figure('Position',[500 100 1000 700])
    subplot(2,1,1);
    plot(P.buf:P.dt:P.tmax-P.buf+180*P.dt,trimdsignals.X+max(trimdsignals.Y)*1.2,'b',...
         P.buf:P.dt:P.tmax-P.buf+180*P.dt,trimdsignals.Y,'r','linewidth',2)
    legend('trimmed signal 1','trimmed signal 2')
    xlim([P.buf P.tmax-P.buf+0.3])
    subplot(2,1,2);
    plot(lags,corrCo)
    title('Cross coorelation')
    saveas(gcf,'crossCo','png');
    
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

function [C,lags,trimdsignals] = CrossThem(x1,x2,P)

% trim the signals to only the important stuff
X = x1(P.buf/P.dt:end-P.buf/P.dt);
Y = x2(P.buf/P.dt:end-P.buf/P.dt);

% do the xcorr
[C,lags] = xcorr(X,Y,'normalized');

% adjust scale of time step
lags = lags*P.dt;

% prepare for output
trimdsignals.X = X;
trimdsignals.Y = Y;

% debuggin figure
if 0 
    figure
    plot(X); hold on
    plot(Y)
end

end
