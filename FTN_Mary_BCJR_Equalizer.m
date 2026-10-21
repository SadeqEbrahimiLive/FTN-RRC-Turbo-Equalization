function L_e_bits_int = FTN_Mary_BCJR_Equalizer( ...
    y, h, noiseVarComplex, mapping, L_a_bits_int, iter)
% FTN_MARY_BCJR_EQUALIZER
%   Log-domain SISO BCJR equalizer for generic M-ary FTN ISI channels.
%
% LLR convention used everywhere:
%   L = log(P(bit=0)/P(bit=1)).
%
% INPUTS:
%   y               : [Ns x 1] received samples
%   h               : [Lch+1 x 1] causal detector-channel taps
%   noiseVarComplex : E{|n|^2} for circular complex Gaussian noise
%   mapping         : constellation/bit-label structure
%   L_a_bits_int    : [Ns x k] a-priori coded-bit LLRs, interleaved order
%   iter            : iteration index (kept for interface compatibility)
%
% OUTPUT:
%   L_e_bits_int    : [Ns x k] EXTRINSIC coded-bit LLRs
%
% The initial ISI state is treated as unknown/equiprobable because the main
% simulation uses a random warm-up sequence before the information block.

if size(y,2) ~= 1
    error('y must be a column vector.');
end
h = h(:);

Ns  = length(y);
Lch = length(h)-1;
M   = mapping.M;
k   = mapping.bitsPerSym;

if size(L_a_bits_int,1) ~= Ns || size(L_a_bits_int,2) ~= k
    error('L_a_bits_int must be Ns-by-k.');
end
if ~(isscalar(noiseVarComplex) && noiseVarComplex > 0)
    error('noiseVarComplex must be a positive scalar E{|n|^2}.');
end

constellation = mapping.constellation(:).';
symBits       = mapping.bits;
numStates     = M^max(Lch,0);

if numStates > 2^12
    warning('BCJR has %d states; runtime may be very large.', numStates);
end

%% ===================== 1. State representation =========================
% State = [a_{n-1}, a_{n-2}, ..., a_{n-Lch}]
if Lch > 0
    stateSymIdx = zeros(numStates,Lch);
    for s = 0:numStates-1
        tmp = s;
        for ell = Lch:-1:1
            digit = mod(tmp,M);
            tmp = floor(tmp/M);
            stateSymIdx(s+1,ell) = digit+1;
        end
    end
else
    stateSymIdx = zeros(1,0);
    numStates = 1;
end

%% ===================== 2. Trellis branches =============================
nextState  = zeros(numStates,M);
outSymHist = cell(numStates,M);
yHatBranch = zeros(numStates,M);

for sPrev = 1:numStates
    past = stateSymIdx(sPrev,:);
    for m = 1:M
        if Lch > 0
            newHist = [m, past(1:end-1)];
            sNext0 = 0;
            for ell = 1:Lch
                sNext0 = sNext0*M + (newHist(ell)-1);
            end
            sNext = sNext0+1;
        else
            newHist = [];
            sNext = 1;
        end

        nextState(sPrev,m) = sNext;
        histIdx = [m,past];
        outSymHist{sPrev,m} = histIdx;
        aVec = constellation(histIdx);
        yHatBranch(sPrev,m) = sum(h(:).'.*aVec);
    end
end

% Incoming branches for each next state
inPrevState = cell(numStates,1);
inInputSym  = cell(numStates,1);
for sPrev = 1:numStates
    for m = 1:M
        sNext = nextState(sPrev,m);
        inPrevState{sNext}(end+1) = sPrev; 
        inInputSym{sNext}(end+1)  = m;    
    end
end

%% ===================== 3. Symbol priors ================================
% L=log(P0/P1), therefore log P(b) = 0.5*L*(1-2b) + constant.
logP_sym = zeros(Ns,M);
for n = 1:Ns
    Li = L_a_bits_int(n,:);
    for m = 1:M
        b = symBits(m,:);
        logP_sym(n,m) = 0.5*sum(Li.*(1-2*b));
    end
end

%% ===================== 4. BCJR alpha/beta ==============================
logAlpha = -inf(Ns+1,numStates);
logBeta  = -inf(Ns+1,numStates);

% Unknown warm-up symbols -> equiprobable initial detector state.
logAlpha(1,:) = 0;
% Unterminated end -> all terminal states equiprobable.
logBeta(Ns+1,:) = 0;

invNoiseVar = 1/noiseVarComplex;

% Forward recursion
for n = 1:Ns
    branchMetric = -abs(y(n)-yHatBranch).^2*invNoiseVar + ...
                   repmat(logP_sym(n,:),numStates,1);

    for sNext = 1:numStates
        prevStates = inPrevState{sNext};
        inSyms = inInputSym{sNext};
        vals = -inf(1,numel(prevStates));
        for q = 1:numel(prevStates)
            sPrev = prevStates(q);
            m = inSyms(q);
            vals(q) = logAlpha(n,sPrev) + branchMetric(sPrev,m);
        end
        logAlpha(n+1,sNext) = logsumexp_vec(vals);
    end

    c = max(logAlpha(n+1,:));
    if isfinite(c)
        logAlpha(n+1,:) = logAlpha(n+1,:)-c;
    end
end

% Backward recursion
for n = Ns:-1:1
    branchMetric = -abs(y(n)-yHatBranch).^2*invNoiseVar + ...
                   repmat(logP_sym(n,:),numStates,1);

    for sPrev = 1:numStates
        vals = -inf(1,M);
        for m = 1:M
            sNext = nextState(sPrev,m);
            vals(m) = branchMetric(sPrev,m) + logBeta(n+1,sNext);
        end
        logBeta(n,sPrev) = logsumexp_vec(vals);
    end

    c = max(logBeta(n,:));
    if isfinite(c)
        logBeta(n,:) = logBeta(n,:)-c;
    end
end

%% ===================== 5. Bit posterior/extrinsic LLRs =================
L_post = zeros(Ns,k);
Nbranch = numStates*M;

for n = 1:Ns
    branchMetric = -abs(y(n)-yHatBranch).^2*invNoiseVar + ...
                   repmat(logP_sym(n,:),numStates,1);

    % Compute the complete branch APP metric once for this time index.
    branchAPP = -inf(numStates,M);
    for sPrev = 1:numStates
        for m = 1:M
            sNext = nextState(sPrev,m);
            branchAPP(sPrev,m) = logAlpha(n,sPrev) + ...
                                 branchMetric(sPrev,m) + ...
                                 logBeta(n+1,sNext);
        end
    end

    for ib = 1:k
        vals0 = -inf(Nbranch,1);
        vals1 = -inf(Nbranch,1);
        n0 = 0;
        n1 = 0;

        for sPrev = 1:numStates
            for m = 1:M
                if symBits(m,ib)==0
                    n0=n0+1;
                    vals0(n0)=branchAPP(sPrev,m);
                else
                    n1=n1+1;
                    vals1(n1)=branchAPP(sPrev,m);
                end
            end
        end

        L_post(n,ib) = logsumexp_vec(vals0(1:n0)) - ...
                       logsumexp_vec(vals1(1:n1));
    end
end

% Equalizer extrinsic = posterior - equalizer a-priori.
L_e_bits_int = L_post - L_a_bits_int;

% Numerical clipping only; no iteration-dependent artificial scaling.
LLRmax = 50;
L_e_bits_int = min(max(L_e_bits_int,-LLRmax),LLRmax);

% Interface compatibility.
if nargin < 6
    iter = 1;
end

end

%% ======================================================================
function s = logsumexp_vec(x)
if isempty(x)
    s = -inf;
    return;
end
m = max(x);
if isinf(m)
    s = m;
else
    s = m + log(sum(exp(x-m)));
end
end
