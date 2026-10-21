# FTN-RRC Turbo Equalization (MATLAB)

A research-oriented MATLAB simulation of **faster-than-Nyquist (FTN) signaling** with **root-raised-cosine (RRC / SRRC) pulse shaping**, M-ary QAM/PSK modulation, a **reduced-memory soft-input/soft-output (SISO) BCJR equalizer**, and iterative exchange of extrinsic information with a rate-1/2 convolutional-code BCJR decoder.

The project contains **RRC/SRRC only**; it does not implement CES pulse shaping.

> **Research status:** This is a simulation implementation, not an experimentally validated modem. Review the system model and run your own validation before using numerical results in a paper. A MATLAB runtime was not available to test the packaged files.

## Quick start

**Requirements**

- MATLAB with support for local functions in scripts (R2016b or later).
- **Communications Toolbox** for `rcosdesign`, `poly2trellis`, `qammod` / `pskmod`, and `de2bi`.
- **Signal Processing Toolbox** for `freqz`.
- Enough memory and CPU time for the BCJR equalizer (see [Performance](#performance)).

**Run**

1. Download/clone the repository and keep the three `.m` files in the **same folder**.
2. Open that folder in MATLAB (or set it as the Current Folder).
3. Run the main script from the Command Window:

   ```matlab
   FTN_Mary_Turbo_RRC_main
   ```

   Alternatively, from a terminal with MATLAB on your PATH:

   ```bash
   matlab -batch "FTN_Mary_Turbo_RRC_main"
   ```

4. Inspect the Command Window for ISI and spectral metrics, and the two generated figures:
   - **Normalized pulse PSD vs. FTN folded spectrum**.
   - **Post-decoder BER vs. information-bit `Eb/N0`**, after the configured final turbo iteration.

The BER results are also available in the MATLAB workspace as `BER_turbo` (rows: SNR points; columns: turbo iterations), `EbN0dB_vec`, `NbitsBER`, and `NerrBER`.

## Files

| File | Purpose |
| --- | --- |
| [`FTN_Mary_Turbo_RRC_main.m`](FTN_Mary_Turbo_RRC_main.m) | Set up the simulation, generate the RRC pulse, analyze FTN interference/spectrum, form an equivalent Forney-factor channel, run Monte Carlo BER estimation, and plot results. **Run this file.** |
| [`FTN_Mary_BCJR_Equalizer.m`](FTN_Mary_BCJR_Equalizer.m) | Log-domain SISO M-ary BCJR equalizer for the truncated ISI channel. Outputs *extrinsic* coded-bit log-likelihood ratios (LLRs). |
| [`Siso_Conv_Decoder_r12.m`](Siso_Conv_Decoder_r12.m) | Log-domain BCJR decoder for the binary rate-1/2 convolutional code. Outputs systematic/parity *extrinsic* LLRs and information-bit a-posteriori LLRs. |
| [`CITATION.cff`](CITATION.cff) | Machine-readable software citation metadata used by GitHub's **Cite this repository** feature. |
| [`LICENSE`](LICENSE) | Proposed MIT license; **replace the copyright-owner placeholder and confirm you want this license before publishing**. |

## Default experiment

Edit parameters in **Section 1** of `FTN_Mary_Turbo_RRC_main.m`.

| Parameter | Default | Meaning |
| --- | --- | --- |
| `modType` | `'qam'` | `'qam'` or `'psk'` |
| `M` | `8` | Modulation alphabet size; use a supported power of two so `k = log2(M)` is an integer |
| `Rc` | `1/2` | Code rate (the implemented encoder/decoder are fixed at 1/2) |
| `Rs` | `1e6` | Nominal symbol rate in symbols/s |
| `SPS` | `10` | Samples per nominal symbol |
| `tau` | `0.7` | Requested FTN time-packing factor; `tau < 1` indicates faster-than-Nyquist packing |
| `beta` | `0.15` | RRC roll-off factor |
| `spanSRRC` | `16` | Total RRC filter span in nominal symbol periods |
| `Lch` | `3` | Equalizer channel memory in symbols (`M^Lch` BCJR states) |
| `Imax` | `5` | Number of equalizer/decoder turbo iterations |
| `EbN0dB_vec` | `0:2:14` | SNR sweep in dB |
| `infoBitsPerFrame` | `3000` | Information bits in each Monte Carlo frame |
| `minErrsFinal` | `200` | Target errors **at the final iteration** before stopping each SNR point |
| `maxInfoBits` | `1e6` | Maximum information bits per SNR point |
| `feedbackDamping` | `1.0` | Scale applied to decoder-to-equalizer extrinsic feedback |

Because FTN samples are placed on a discrete grid, the effective packing factor is

```matlab
L_FTN = round(tau*SPS);
tauEff = L_FTN/SPS;
```

`tauEff`, not necessarily the requested `tau`, is what the code uses. Choose `SPS` and `tau` accordingly. Keep `spanSRRC * SPS` even, as required by `rcosdesign`.

### Short smoke test

The default sweep can be very slow. For a **small functional check only**, temporarily change these values in the main script:

```matlab
EbN0dB_vec     = 4;  % one SNR point
Imax           = 2;  % two iterations
Lch            = 2;  % 8^2 = 64 BCJR states, not the default 512
infoBitsPerFrame = 300;
minErrsFinal   = 20;
maxInfoBits    = 300;
```

**Do not use smoke-test results for publication:** they are based on too few bits and a different detector memory. Restore the defaults (or specify and justify a proper experiment) for research-grade simulations.

### Other experiments

```matlab
modType = 'psk';  M = 8;   % 8-PSK (keep M a power of two)
beta = 0.25;              % a different RRC roll-off
EbN0dB_vec = 0:1:12;     % finer SNR spacing
Lch = 2;                 % faster, lower-complexity detector
```

Change **one design choice at a time** when comparing BER curves; save the associated parameters and random seed for each reported result.

## Processing chain

```text
Random information bits
    -> rate-1/2 systematic convolutional encoding (generators [4 7], octal)
    -> random bit interleaver
    -> Gray-labeled M-QAM or M-PSK mapping
    -> equivalent causal FTN/Forney-factor ISI channel + complex AWGN
    -> reduced-memory M-ary SISO BCJR equalizer
    -> deinterleaver
    -> systematic/parity SISO convolutional BCJR decoder
    -> interleaved extrinsic feedback (repeat Imax times)
    -> final information-bit decisions and BER
```

**RRC and channel model.** The script creates a unit-energy RRC pulse with `rcosdesign(...,'sqrt')`. It samples the pulse autocorrelation at FTN spacing and computes a minimum-phase spectral factor. The BER simulation passes mapped symbols through the resulting **equivalent causal discrete-time spectral-factor channel**; it does *not* explicitly upsample, pulse-shape, matched-filter, and sample the waveform at every Monte Carlo frame. The full retained factor is used to generate observations, whereas the BCJR uses only `Lch + 1` taps. Energy in the omitted tail is approximated as extra branch uncertainty:

```matlab
noiseVarEff = N0 + EsTx*residualEnergy;
```

This is an approximation; the omitted ISI is not actually independent Gaussian noise. The spectral-factor truncation and detector memory therefore matter when interpreting results.

**Code/decoder.** The implemented feed-forward binary code uses `poly2trellis(3,[4 7])` with a systematic first output and a parity second output. The parity is generated from the current and two preceding input bits. Frames start in the zero state, but **no termination bits** are appended, so the decoder allows any final state.

**LLR convention.** Both helper functions use `L = log(P(bit=0)/P(bit=1))`; the final hard decision is `L < 0`. Each module passes **extrinsic** information to the other module. The decoder returns extrinsic feedback for **both** systematic and parity coded bits.

**SNR definition.** The script uses `EbInfo = EsTx / (Rc*k)` and complex noise variance `N0 = EbInfo / 10^(EbN0dB/10)`. `Eb/N0` is referenced to information bits, not coded bits.

## Diagnostics and interpreting results

The command output includes:

- FTN-spaced autocorrelation taps and main-to-ISI power ratio.
- `B99`: numerically estimated **full two-sided 99%-occupied bandwidth** of the pulse response.
- Gross/net spectral efficiency based on nominal/FTN bit rates and `B99`.
- Out-of-cell alias energy ratio (OC-AER), folded-spectrum ripple, and spectral condition metrics.
- Spectral-factor energy capture, truncated-detector taps, and tail-energy percentage.
- Bits simulated and errors after **each turbo iteration** at every SNR point.

For SNR points with **zero observed final-iteration errors**, the BER plot displays the approximate 95% upper confidence bound `3/N`, shown with a downward triangle. **That marker is not a measured nonzero BER.** The default Monte Carlo loop stops when the final-iteration errors reach `minErrsFinal` or `maxInfoBits` is reached. It does not guarantee comparable precision at every SNR point.

## Performance and limitations

- Detector state count grows as **`M^Lch`**. With the defaults, `8^3 = 512` states; each state has eight outgoing branches. Multiple iterations and large frames can make this MATLAB implementation **computationally expensive**.
- This is a **reduced-memory** detector and uses an effective-noise approximation for residual ISI; it is not full-memory optimal BCJR of the original pulse-shaped physical waveform.
- The simulation uses AWGN and an equivalent ISI factor channel; it does not model carrier offset, timing recovery, frequency-selective fading, RF hardware, or synchronization.
- Changing the convolutional code in `trellis` alone is **not** sufficient: the handwritten encoder and the decoder must stay consistent. `Rc = 1/2` is currently fixed by the implementation.
- The project does not include automated MATLAB tests or published reference BER datasets. Before publishing numerical claims, verify toolbox availability, reproduce runs, document MATLAB version/hardware, and perform appropriate Monte Carlo convergence checks.
- MATLAB/Communications Toolbox/Signal Processing Toolbox are MathWorks products and are not included in this repository.

## Citation

**If you use or extend this software in a paper, please cite the repository or, preferably, its DOI-archived release.** The root-level [`CITATION.cff`](CITATION.cff) is configured to make GitHub show **Cite this repository** after it is pushed.

## License

A proposed **MIT** license is included for broad software reuse. **Before publishing, change `Sadeq Ebrahimi` in `LICENSE` to the true copyright holder(s), confirm that you have rights to license all included code, and decide whether MIT fits your intentions.** MIT asks redistributors to preserve the license/copyright notice; it does **not** make academic citation a legal condition. `CITATION.cff` and the request above facilitate scholarly attribution.

## Feedback / contributions

Bug reports, reproducible performance comparisons, documentation improvements, and pull requests are welcome once the repository is published. For numerical issues, include MATLAB release, toolbox versions, parameter values, RNG configuration, and console output.
