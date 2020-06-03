%% Testing cross coorelation with a toy model
% This simulation aims to test out predictions about using cross
% coorelation on two possibly related neuronal recordings

% First the code generates fictive spike trains. Either poisson proccesses
% or repetative regular spikes with jitter. The user can play with all the
% parameters in the %% Parameters section
clear; close all;

%% Parameters

P.tmax = 20;   % maximum simulation time in secconds
P.dt   = 0.001; % time step

P.traintoggle = 2; % the switch for spike timings, 
                     % If 1, the spike trains will be regular spaced inputs 
                     % with some jitter.
                     % If 2, one spikes train will be a poisson process and
                     % the other will be a copy, with some jitter
                     % If 3, both trains will be different poissons

% regular spikes settings                                           
P.N    = 10;   % number of neurons
P.p    = 1;    % period of spiking
P.vp   = 2;    % STD of spiking

% poisson proccess settings
P.lambda = 0.2; % average timing between spikes in seconds 
P.jit    = 0.1; % std to jitter spike timings

P.buf  = 5;    % amount of buffer time before/after spiking. Recommend a 
               % healthy window

%% Presim calculations and initializations

P.t    = 0:P.dt:P.tmax-P.dt; % time array
P.tct  = [0:P.dt:P.tmax-P.dt+(180*P.dt)]'; % time array for convoluted 
                                           % signal (it's longer due to 
                                           % the convolution).

P.nsteps  = P.tmax/P.dt;
P.nspikes = P.tmax/P.p - (P.buf/P.p-1) - (P.buf/P.p);

if  2*P.vp>=P.buf
    warning(['A period variance (P.vp) larger than or close too the front buffer ' ...
             '(P.buf) may push a spike before the simulation start, causing an error. '...
             ' Recommend P.vp<<P.buf'])
end

%% Blocks

% generate both spike trains
switch P.traintoggle
    case 1
    % if the user wants regularly spaced spikes
    [sumt1,tra1] = GenSpikTrain(P);
    [sumt2,tra2] = GenSpikTrain(P);
    case 2
    % if the user wants a poisson distanced spikes
    [sumt1] = GenSpikPoisT(P);
    [sumt2] = JitterTrain(sumt1,P);
    case 3
    [sumt1] = GenSpikPoisT(P);
    [sumt2] = GenSpikPoisT(P);
end

% it's possible t1 and t2 are different sizes if its the poisson trains
% just need to fix that
if length(sumt1)~=length(sumt2) 
    m = min([length(sumt1) length(sumt2)]);
    sumt1 = sumt1(1:m);
    sumt2 = sumt2(1:m);
end

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
    xlabel('Time (s)')
    saveas(gcf,'introduction','png'); %close
    
    
    figure('Position',[500 100 1000 700])
    subplot(2,1,1);
    plot(P.buf:P.dt:P.tmax-P.buf+180*P.dt,trimdsignals.X+max(trimdsignals.Y)*1.2,'b',...
         P.buf:P.dt:P.tmax-P.buf+180*P.dt,trimdsignals.Y,'r','linewidth',2)
    legend('trimmed signal 1','trimmed signal 2')
    xlim([P.buf P.tmax-P.buf+0.3])
    xlabel('Time (s)')
    subplot(2,1,2);
    plot(lags,corrCo)
    title('Cross coorelation')
    xlabel('Lag (s)')
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

function [train] = GenSpikPoisT(P)

    % adjusted lambda, adjusted for the dt and index stuff
    adj_lam = P.lambda/P.dt;

    % initialize train
    train = zeros(length(P.t),1);

    for k = 1:P.tmax/P.dt; if rand<1/adj_lam; train(k) = 1; end; end
    
    % debugging figure
    if 0
        figure
        plot(P.t,train)
        % average distance between spikes = 
        adbs = mean(diff(find(train)))*P.dt
    end
    
end

function [t2] = JitterTrain(t1,P)

    x = find(t1);
    r = round( rand(length(x),1) * ( P.jit/P.dt ) );
    y = x+r;
    
    t2 = zeros(length(t1),1);
    t2(y) = 1;
    
    if 0
        figure
        plot(P.t,t1); hold on;
        plot(P.t,t2);
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
