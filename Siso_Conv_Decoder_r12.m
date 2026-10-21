function [L_sys_ext, L_par_ext, L_sys_post] = Siso_Conv_Decoder_r12( ...
    L_ch_sys, L_ch_par, L_a_sys, trellis, iter)
% SISO_CONV_DECODER_R12
%   Log-domain BCJR decoder for the systematic rate-1/2 convolutional code.
%
% LLR convention used everywhere:
%   L = log(P(bit=0)/P(bit=1)).
%
% INPUTS:
%   L_ch_sys : [K x 1] equalizer EXTRINSIC LLRs for systematic bits
%   L_ch_par : [K x 1] equalizer EXTRINSIC LLRs for parity bits
%   L_a_sys  : [K x 1] optional extra a-priori LLRs on information bits
%              (set to zero in the corrected turbo loop)
%   trellis  : MATLAB convolutional-code trellis
%   iter     : iteration index (kept for interface compatibility)
%
% OUTPUTS:
%   L_sys_ext  : decoder EXTRINSIC LLRs for systematic coded bits
%   L_par_ext  : decoder EXTRINSIC LLRs for parity coded bits
%   L_sys_post : a-posteriori information-bit LLRs for hard decisions
%
% Returning parity extrinsic information is important for bit-interleaved
% M-ary turbo equalization: all coded-bit positions are fed back to the
% M-ary equalizer, not only the systematic positions.

%% ===================== 1. Basic checks =================================
if size(L_ch_sys,2)~=1 || size(L_ch_par,2)~=1 || size(L_a_sys,2)~=1
    error('All LLR inputs must be column vectors.');
end
K = length(L_ch_sys);
if length(L_ch_par)~=K || length(L_a_sys)~=K
    error('LLR vectors must have equal length.');
end

numStates = trellis.numStates;
if trellis.numInputSymbols~=2 || trellis.numOutputSymbols~=4
    error('This decoder assumes a binary-input rate-1/2 trellis.');
end

nextStates = trellis.nextStates + 1;
outputs = trellis.outputs;

%% ===================== 2. Trellis branch labels ========================
outNext = zeros(numStates,2);
outU    = zeros(numStates,2);
outSys  = zeros(numStates,2);
outPar  = zeros(numStates,2);

for s = 1:numStates
    for uIdx = 1:2
        outNext(s,uIdx) = nextStates(s,uIdx);
        outU(s,uIdx) = uIdx-1;
        ob = de2bi(outputs(s,uIdx),2,'left-msb');
        outSys(s,uIdx) = ob(1);
        outPar(s,uIdx) = ob(2);
    end
end

inPrev = cell(numStates,1);
inU    = cell(numStates,1);
inSys  = cell(numStates,1);
inPar  = cell(numStates,1);

for sPrev = 1:numStates
    for uIdx = 1:2
        sNext = outNext(sPrev,uIdx);
        inPrev{sNext}(end+1)=sPrev; 
        inU{sNext}(end+1)=outU(sPrev,uIdx);
        inSys{sNext}(end+1)=outSys(sPrev,uIdx);
        inPar{sNext}(end+1)=outPar(sPrev,uIdx);
    end
end

L_ch_sys = L_ch_sys(:).';
L_ch_par = L_ch_par(:).';
L_a_sys  = L_a_sys(:).';

%% ===================== 3. Alpha/beta recursions ========================
logAlpha = -inf(K+1,numStates);
logBeta  = -inf(K+1,numStates);

% Encoder starts in the all-zero state.
logAlpha(1,1)=0;
% No tail bits are appended, so final state is unknown.
logBeta(K+1,:)=0;

for kk = 1:K
    for sNext = 1:numStates
        prev = inPrev{sNext};
        vals = -inf(1,numel(prev));
        for q=1:numel(prev)
            sPrev = prev(q);
            u  = inU{sNext}(q);
            bs = inSys{sNext}(q);
            bp = inPar{sNext}(q);

            gamma = 0.5*L_ch_sys(kk)*(1-2*bs) + ...
                    0.5*L_ch_par(kk)*(1-2*bp) + ...
                    0.5*L_a_sys(kk)*(1-2*u);

            vals(q)=logAlpha(kk,sPrev)+gamma;
        end
        logAlpha(kk+1,sNext)=logsumexp_vec(vals);
    end

    c=max(logAlpha(kk+1,:));
    if isfinite(c)
        logAlpha(kk+1,:)=logAlpha(kk+1,:)-c;
    end
end

for kk=K:-1:1
    for sPrev=1:numStates
        vals=-inf(1,2);
        for uIdx=1:2
            sNext=outNext(sPrev,uIdx);
            u = outU(sPrev,uIdx);
            bs=outSys(sPrev,uIdx);
            bp=outPar(sPrev,uIdx);

            gamma = 0.5*L_ch_sys(kk)*(1-2*bs) + ...
                    0.5*L_ch_par(kk)*(1-2*bp) + ...
                    0.5*L_a_sys(kk)*(1-2*u);

            vals(uIdx)=gamma+logBeta(kk+1,sNext);
        end
        logBeta(kk,sPrev)=logsumexp_vec(vals);
    end

    c=max(logBeta(kk,:));
    if isfinite(c)
        logBeta(kk,:)=logBeta(kk,:)-c;
    end
end

%% ===================== 4. Systematic and parity APP LLRs ===============
L_sys_post=zeros(K,1);
L_par_post=zeros(K,1);

for kk=1:K
    sys0=[]; sys1=[];
    par0=[]; par1=[];

    for sPrev=1:numStates
        for uIdx=1:2
            sNext=outNext(sPrev,uIdx);
            u = outU(sPrev,uIdx);
            bs=outSys(sPrev,uIdx);
            bp=outPar(sPrev,uIdx);

            gamma = 0.5*L_ch_sys(kk)*(1-2*bs) + ...
                    0.5*L_ch_par(kk)*(1-2*bp) + ...
                    0.5*L_a_sys(kk)*(1-2*u);

            val=logAlpha(kk,sPrev)+gamma+logBeta(kk+1,sNext);

            if u==0
                sys0(end+1)=val;
            else
                sys1(end+1)=val;
            end

            if bp==0
                par0(end+1)=val;
            else
                par1(end+1)=val;
            end
        end
    end

    L_sys_post(kk)=logsumexp_vec(sys0)-logsumexp_vec(sys1);
    L_par_post(kk)=logsumexp_vec(par0)-logsumexp_vec(par1);
end

%% ===================== 5. Decoder extrinsic information =================
% The equalizer outputs are already extrinsic, so the main turbo loop uses
% L_a_sys=0.  Remove each coded bit's own input before feeding back.
L_sys_ext=L_sys_post-L_ch_sys(:)-L_a_sys(:);
L_par_ext=L_par_post-L_ch_par(:);

LLRmax=50;
L_sys_ext=min(max(L_sys_ext,-LLRmax),LLRmax);
L_par_ext=min(max(L_par_ext,-LLRmax),LLRmax);
L_sys_post=min(max(L_sys_post,-LLRmax),LLRmax);

if nargin<5
    iter=1;
end

end

%% ======================================================================
function s=logsumexp_vec(x)
if isempty(x)
    s=-inf;
    return;
end
m=max(x);
if isinf(m)
    s=m;
else
    s=m+log(sum(exp(x-m)));
end
end
