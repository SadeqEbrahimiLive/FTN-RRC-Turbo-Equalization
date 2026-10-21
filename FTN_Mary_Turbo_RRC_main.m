%% FTN_Mary_Turbo_RRC_main.m
% FTN turbo equalization with square-root raised-cosine (SRRC/RRC) pulses.
% Supports M-QAM / M-PSK, FTN ISI and spectrum diagnostics, Forney
% spectral factorization, reduced-memory BCJR, turbo decoding, and BER.
%
% Requires:
%   - FTN_Mary_BCJR_Equalizer.m
%   - Siso_Conv_Decoder_r12.m
%   - Communications Toolbox (rcosdesign, poly2trellis, qammod/pskmod, de2bi)

clear; clc; close all;

%% 1. System parameters
modType = 'qam';
M       = 8;
k       = log2(M);
Rc      = 1/2;

Rs      = 1e6;              % nominal symbol rate [symbols/s]
RbGross = Rs*k;              % coded-bit rate for tau = 1
RbInfo  = Rc*RbGross;        % information-bit rate for tau = 1

SPS     = 10;               % samples per nominal symbol
fs      = SPS*Rs;           % sample frequency [Hz]
tau     = 0.7;              % FTN packing factor
L_FTN   = round(tau*SPS);   % integer sample spacing
tauEff  = L_FTN/SPS;

% Square-root raised-cosine pulse parameters.
spanSRRC = 16;              % filter span in nominal symbols
beta     = 0.15;            % roll-off factor

% Turbo code and receiver parameters.
trellis = poly2trellis(3,[4 7]);
Imax    = 5;
Lch     = 3;               % BCJR memory: M^Lch states

EbN0dB_vec = 0:2:14;
nSNR       = numel(EbN0dB_vec);

% Monte-Carlo stopping criteria.
infoBitsPerFrame = 3000;
minErrsFinal     = 200;
maxInfoBits      = 1e6;

if mod(2*infoBitsPerFrame,k) ~= 0
    error('2*infoBitsPerFrame must be divisible by log2(M).');
end

% Spectral factorization parameters.
NfftFactor       = 16384;
spectralFloorRel = 1e-12;
factorEnergyKeep = 1-1e-8;
maxFactorTaps    = 128;

% Spectral diagnostic parameters.
metricMask  = 0.99;
NfftMetrics = 16384;

feedbackDamping = 1.0;      % decoder-to-equalizer feedback damping

fprintf('Pulse: SRRC, beta = %.3f, span = %d symbols\n',beta,spanSRRC);
fprintf('Rs = %.3e sym/s, gross Rb = %.3e bit/s, net Rb = %.3e bit/s\n', ...
    Rs,RbGross,RbInfo);
fprintf('fs = %.3e Hz, SPS = %d, tau = %.3f, L_FTN = %d samples\n', ...
    fs,SPS,tauEff,L_FTN);
fprintf('Modulation: %s, M=%d, k=%d, Rc=%.3f\n',upper(modType),M,k,Rc);
fprintf('BCJR detector memory Lch=%d (%d states)\n',Lch,M^Lch);

mapping = create_mapping_struct(modType,M);

%% 2. SRRC pulse design
h_tx = rcosdesign(beta,spanSRRC,SPS,'sqrt').';
h_tx = h_tx(:);
h_tx = h_tx/sqrt(sum(abs(h_tx).^2));

fprintf('Pulse length=%d, energy=%.12f\n', ...
    numel(h_tx),sum(abs(h_tx).^2));

%% 3. FTN pulse/channel analysis
% 3.1 Sampled autocorrelation and ISI
Kmax = floor((numel(h_tx)-1)/L_FTN);
[gpos,mainPow,isiPow,mainIsi_dB] = sampled_autocorr(h_tx,L_FTN,Kmax);

fprintf('Autocorrelation taps g[k], k=0..min(7,Kmax):\n');
disp(gpos(1:min(8,numel(gpos))).');
fprintf('Main power = %.6g, off-center ISI power = %.6g, main/ISI = %.2f dB\n', ...
    mainPow,isiPow,mainIsi_dB);

% 3.2 Spectrum, 99%% occupied bandwidth, and spectral efficiency
met = evaluate_pulse_metrics(h_tx,fs,RbGross,RbInfo,tauEff, ...
    L_FTN,NfftMetrics,metricMask);

fprintf('B99(full)        = %.6f MHz\n',met.BW99/1e6);
fprintf('SE gross Nyq/FTN = %.4f / %.4f bit/s/Hz\n', ...
    RbGross/met.BW99,met.SE_gross);
fprintf('SE net   Nyq/FTN = %.4f / %.4f bit/s/Hz\n', ...
    RbInfo/met.BW99,met.SE_net);
fprintf('OC-AER           = %.3f dB\n',met.OCAER_dB);
fprintf('Fold PA ripple   = %.3f dB\n',met.RipplePA_dB);
fprintf('Arel=min/mean    = %.6e\n',met.Arel);
fprintf('Brel=max/mean    = %.6e\n',met.Brel);
fprintf('kappa_inf/spec   = %.6g (%.3f dB)\n', ...
    met.kappa_spec,met.DynRange_dB);

% 3.3 Minimum-phase Forney factor / reduced-memory BCJR model
[fFactor,~,SrawMin,SrawMax] = minphase_factor_from_gpos( ...
    gpos,NfftFactor,spectralFloorRel);
[h_full,capturedFull] = truncate_factor_energy( ...
    fFactor,factorEnergyKeep,maxFactorTaps);

h_det = h_full(1:min(Lch+1,numel(h_full)));
if numel(h_det) < Lch+1
    h_det(end+1:Lch+1) = 0;
end

captureDet       = sum(abs(h_det).^2)/sum(abs(h_full).^2);
residualDet      = max(0,1-captureDet);
residualEnergy   = max(0,sum(abs(h_full).^2)-sum(abs(h_det).^2));
kappaFactor      = SrawMax/max(SrawMin,realmin);
kappaMismatchPct = 100*abs(kappaFactor-met.kappa_spec) / ...
    max(met.kappa_spec,realmin);

fprintf('Forney factor: Smin=%.4e, Smax=%.4e, max/min=%.6g (%.3f dB)\n', ...
    SrawMin,SrawMax,kappaFactor,10*log10(kappaFactor));
fprintf('Direct-fold vs Forney kappa difference = %.3f%%\n',kappaMismatchPct);
fprintf('Factor taps=%d (captured %.8f of factor energy)\n', ...
    numel(h_full),capturedFull);
fprintf('BCJR taps=%d (Lch=%d), captures %.3f%%, residual tail=%.3f%%\n', ...
    numel(h_det),Lch,100*captureDet,100*residualDet);
fprintf('Detector taps h_det:\n');
disp(h_det(:).');

if kappaMismatchPct > 5
    warning('SRRC folded-spectrum and Forney-factor condition ratios differ by %.1f%%.', ...
        kappaMismatchPct);
end
if captureDet < 0.95
    warning('SRRC BCJR memory captures only %.1f%% of factor energy.', ...
        100*captureDet);
end

%% 4. SRRC PSD and FTN folded spectrum
figure('Color','w');
plot(met.f/1e6,10*log10(met.Pow/max(met.Pow)+1e-14), ...
    'LineWidth',1.2);
hold on;
Sf = met.fold_norm/max(met.fold_norm);
plot(met.fold_f/1e6,10*log10(Sf+1e-14),'--','LineWidth',1.5);
grid on;
xlim([-2.5 2.5]);
ylim([-90 3]);
xlabel('Frequency (MHz)');
ylabel('Normalized power (dB)');
title(sprintf('SRRC: PSD and FTN folded spectrum, \\tau=%.2f',tauEff));
legend('Pulse |H(f)|^2','Folded spectrum','Location','best');
hold off;

%% 5. BER Monte Carlo for SRRC only
BER_turbo = nan(nSNR,Imax);
NbitsBER  = zeros(nSNR,1);
NerrBER   = zeros(nSNR,Imax);

for iSNR = 1:nSNR
    EbN0dB  = EbN0dB_vec(iSNR);
    EbN0lin = 10^(EbN0dB/10);
    fprintf('\n[SRRC] Eb/N0 = %.1f dB\n',EbN0dB);

    rng(10000+iSNR,'twister');
    errCount = zeros(Imax,1);
    bitCount = 0;
    frame    = 0;

    while (errCount(end) < minErrsFinal) && (bitCount < maxInfoBits)
        frame = frame+1;

        Kinfo = min(infoBitsPerFrame,maxInfoBits-bitCount);
        while mod(2*Kinfo,k) ~= 0
            Kinfo = Kinfo-1;
        end
        if Kinfo <= 0
            break;
        end
        Ns = 2*Kinfo/k;

        % 5.1 Systematic rate-1/2 convolutional encoder
        u = randi([0 1],Kinfo,1);
        uPad = [0;0;u];
        v = mod(uPad(3:end)+uPad(2:end-1)+uPad(1:end-2),2);

        coded = zeros(2*Kinfo,1);
        coded(1:2:end) = u;
        coded(2:2:end) = v;
        posSys = 1:2:2*Kinfo;
        posPar = 2:2:2*Kinfo;

        % 5.2 Interleaving and M-ary symbol mapping
        piBits = randperm(2*Kinfo);
        piDeint = zeros(1,2*Kinfo);
        piDeint(piBits) = 1:2*Kinfo;

        codedInt = coded(piBits);
        bitsSym = reshape(codedInt,k,Ns).';

        symTx = zeros(Ns,1);
        for n = 1:Ns
            idxSym = find(all(mapping.bits==bitsSym(n,:),2),1,'first');
            symTx(n) = mapping.constellation(idxSym);
        end

        EsTx   = mean(abs(symTx).^2);
        EbInfo = EsTx/(Rc*k);
        N0     = EbInfo/EbN0lin;

        % 5.3 Equivalent causal Forney-factor channel and noise
        Npre = max(20,numel(h_full)+Lch);
        preIdx = randi([1 M],Npre,1);
        preSym = mapping.constellation(preIdx);
        x = [preSym(:);symTx];
        yAll = filter(h_full,1,x);
        yClean = yAll(Npre+1:Npre+Ns);

        noise = sqrt(N0/2)*(randn(Ns,1)+1j*randn(Ns,1));
        y = yClean+noise;

        % Account for factor taps omitted from the BCJR branch model.
        noiseVarEff = N0+EsTx*residualEnergy;

        % 5.4 Iterative BCJR equalization and turbo decoding
        L_a_bits = zeros(2*Kinfo,1);
        for iter = 1:Imax
            L_a_int_vec = L_a_bits(piBits);
            L_a_int = reshape(L_a_int_vec,k,Ns).';

            L_e_eq_int = FTN_Mary_BCJR_Equalizer( ...
                y,h_det,noiseVarEff,mapping,L_a_int,iter);
            L_e_eq_vec = reshape(L_e_eq_int.',[],1);
            L_e_eq = L_e_eq_vec(piDeint);

            L_eq_sys = L_e_eq(posSys);
            L_eq_par = L_e_eq(posPar);
            L_a_dec = zeros(Kinfo,1);

            [L_dec_ext_sys,L_dec_ext_par,L_dec_post] = ...
                Siso_Conv_Decoder_r12( ...
                    L_eq_sys,L_eq_par,L_a_dec,trellis,iter);

            L_a_bits_new = zeros(2*Kinfo,1);
            L_a_bits_new(posSys) = feedbackDamping*L_dec_ext_sys;
            L_a_bits_new(posPar) = feedbackDamping*L_dec_ext_par;
            L_a_bits = L_a_bits_new;

            uHat = double(L_dec_post<0); % LLR = log(P(bit=0)/P(bit=1))
            errCount(iter) = errCount(iter)+sum(uHat~=u);
        end

        bitCount = bitCount+Kinfo;
        if frame==1 || mod(frame,10)==0
            fprintf('  frame=%d, bits=%d, final errors=%d, BER~%.3g\n', ...
                frame,bitCount,errCount(end),errCount(end)/bitCount);
        end
    end

    BER_turbo(iSNR,:) = errCount.'/max(bitCount,1);
    NbitsBER(iSNR)    = bitCount;
    NerrBER(iSNR,:)   = errCount.';

    fprintf('  FINAL: bits=%d, errors by iter=',bitCount);
    fprintf(' %d',errCount);
    fprintf('\n         BER by iter   =');
    fprintf(' %.4g',BER_turbo(iSNR,:));
    fprintf('\n');
    if errCount(end)==0
        fprintf('         Zero errors: 95%% upper bound ~%.3g (3/N)\n', ...
            3/max(bitCount,1));
    end
end

%% 6. BER figure (SRRC only)
bPlot = BER_turbo(:,Imax);
zeroIdx = NerrBER(:,Imax)==0;
bPlot(zeroIdx) = 3./max(NbitsBER(zeroIdx),1); % 95% upper bounds

figure('Color','w');
hLine = semilogy(EbN0dB_vec,bPlot,'-s','LineWidth',1.8, ...
    'MarkerSize',7);
hold on;
if any(zeroIdx)
    semilogy(EbN0dB_vec(zeroIdx),bPlot(zeroIdx),'v', ...
        'LineStyle','none','Color',hLine.Color, ...
        'MarkerSize',8,'LineWidth',1.3,'HandleVisibility','off');
end

ax = gca;
ax.YScale = 'log';
grid on;
box on;
xlim([min(EbN0dB_vec) max(EbN0dB_vec)]);
ylim([1e-6 3e-1]);
yticks([1e-6 1e-5 1e-4 1e-3 1e-2 1e-1]);
ax.YMinorGrid = 'on';
ax.XMinorGrid = 'off';
ax.FontSize = 12;
ax.FontWeight = 'bold';
ax.LineWidth = 1.1;

xlabel('$E_b/N_0$ (dB)','Interpreter','latex','FontSize',14);
ylabel('Post-decoder BER','FontSize',14);
legend(hLine,{'SRRC'},'Location','southwest','FontSize',11);
title(sprintf(['%d-%s, $R_c=1/2$, $\\tau=%.2f$, ', ...
               '$L_{\\mathrm{ch}}=%d$, %d turbo iterations'], ...
               M,upper(modType),tauEff,Lch,Imax), ...
      'Interpreter','latex','FontSize',13);
hold off;

fprintf('\nDone.\n');

%% SUPPORT FUNCTIONS
function mapping = create_mapping_struct(modType,M)
k = log2(M);
symIdx = (0:M-1).';
switch lower(modType)
    case 'qam'
        const = qammod(symIdx,M,'gray','UnitAveragePower',true);
    case 'psk'
        const = pskmod(symIdx,M,0,'gray');
    otherwise
        error('Unsupported modulation type.');
end
mapping.M = M;
mapping.bitsPerSym = k;
mapping.constellation = const(:).';
mapping.bits = de2bi(symIdx,k,'left-msb');
end

function [gpos,mainPow,isiPow,ratio_dB] = sampled_autocorr(h,L,Kmax)
h = h(:);
gpos = zeros(Kmax+1,1);
for kk = 0:Kmax
    sh = kk*L;
    if sh < numel(h)
        gpos(kk+1) = sum(h(1+sh:end).*conj(h(1:end-sh)));
    end
end
mainPow = abs(gpos(1))^2;
isiPow = 2*sum(abs(gpos(2:end)).^2);
ratio_dB = 10*log10(mainPow/max(isiPow,realmin));
end

function met = evaluate_pulse_metrics(p,fs,Rb_gross,Rb_net,rho,L_ftn,nfft,mask)
p = p(:);
p = p/sqrt(sum(abs(p).^2));
[H,w] = freqz(p,1,nfft,'whole');
H = fftshift(H);
w = w-pi;
f = w*(fs/(2*pi));
Pow = abs(H).^2;
df = f(2)-f(1);
Pow = Pow/(sum(Pow)*df+eps);

% Two-sided 99% occupied bandwidth and FTN spectral efficiency.
BW99 = two_sided_bw_cdf(f,Pow,mask);
SE_gross = (Rb_gross/rho)/BW99;
SE_net = (Rb_net/rho)/BW99;

% Out-of-cell alias energy in the FTN Nyquist cell.
f_alias = fs/L_ftn;
fcell = f_alias/2;
Ein = sum(Pow(abs(f)<=fcell))*df;
Eout = sum(Pow(abs(f)>fcell))*df;
OCAER_dB = 10*log10((Eout+eps)/(Ein+eps));

% Folded spectrum with L_ftn unique aliases modulo fs.
Nfold = max(2048,2^nextpow2(4*nfft/L_ftn));
fold_f = linspace(-f_alias/2,f_alias/2,Nfold).';
Sfold = zeros(size(fold_f));
for r = 0:L_ftn-1
    fshift = fold_f+r*f_alias;
    fwrap = wrapToInterval(fshift,-fs/2,fs/2);
    Sfold = Sfold+interp1(f,Pow,fwrap,'linear',0);
end
Sfold = Sfold/L_ftn;

Smean = mean(Sfold);
Smin = min(Sfold);
Smax = max(Sfold);
Arel = Smin/(Smean+eps);
Brel = Smax/(Smean+eps);
if Smin<=0
    kappa_spec = inf;
else
    kappa_spec = Smax/Smin;
end
DynRange_dB = 10*log10(kappa_spec);
RipplePA_dB = 10*log10(Smax/(Smean+eps));

met = struct;
met.BW99 = BW99;
met.SE_gross = SE_gross;
met.SE_net = SE_net;
met.OCAER_dB = OCAER_dB;
met.RipplePA_dB = RipplePA_dB;
met.Arel = Arel;
met.Brel = Brel;
met.kappa_spec = kappa_spec;
met.DynRange_dB = DynRange_dB;
met.fold_f = fold_f;
met.fold_norm = Sfold/(Smean+eps);
met.f = f;
met.Pow = Pow;
end

function BW = two_sided_bw_cdf(f,PSD,mask)
df = f(2)-f(1);
PSD = PSD/(sum(PSD)*df+eps);
cdfv = cumsum(PSD)*df;
low = (1-mask)/2;
high = 1-low;
iL = find(cdfv>=low,1,'first');
iU = find(cdfv>=high,1,'first');
if isempty(iL) || isempty(iU)
    BW = NaN;
else
    BW = abs(f(iU)-f(iL));
end
end

function [fFactor,S,Smin,Smax] = minphase_factor_from_gpos(gpos,Nfft,floorRel)
% Minimum-phase spectral factor of the sampled autocorrelation sequence.
gpos = gpos(:);
K = numel(gpos)-1;
w = 2*pi*(0:Nfft-1).'/Nfft;
S = real(gpos(1))*ones(Nfft,1);
for kk = 1:K
    S = S+2*real(gpos(kk+1)*exp(-1j*w*kk));
end
Smin = min(S);
Smax = max(S);
Sreg = max(S,max(Smax*floorRel,realmin));

c = real(ifft(log(Sreg)));
cmin = zeros(Nfft,1);
cmin(1) = 0.5*c(1);
if mod(Nfft,2)==0
    cmin(2:Nfft/2) = c(2:Nfft/2);
    cmin(Nfft/2+1) = 0.5*c(Nfft/2+1);
else
    cmin(2:(Nfft+1)/2) = c(2:(Nfft+1)/2);
end
F = exp(fft(cmin));
fFactor = ifft(F);
fFactor = fFactor*exp(-1j*angle(fFactor(1)));
end

function [hKeep,captured] = truncate_factor_energy(f,target,maxTaps)
f = f(:);
Etot = sum(abs(f).^2);
cum = cumsum(abs(f).^2)/Etot;
idx = find(cum>=target,1,'first');
if isempty(idx)
    idx = maxTaps;
end
idx = min(idx,maxTaps);
hKeep = f(1:idx);
captured = sum(abs(hKeep).^2)/Etot;
end

function y = wrapToInterval(x,minVal,maxVal)
L = maxVal-minVal;
y = x-L.*floor((x-minVal)./L);
y(y>=maxVal) = minVal;
end
