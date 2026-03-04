# Physiological Signal Preprocessing Guide

## EDA (Electrodermal Activity)

### Noise sources
- **High-frequency noise**: EMI from equipment (computers, lights) — low-amplitude, high-freq oscillations
- **Motion artifacts**: Large, rapid changes from movement/electrode displacement

### Standard preprocessing (from [EDA Guidelines](https://edaguidelines.github.io/analysis/pre), NeuroKit2, pyEDA)

| Step | Method | Parameters |
|------|--------|------------|
| **1. Low-pass filter** | Butterworth 4th order | Cutoff **3–5 Hz** (NeuroKit: 3 Hz, BioSPPy: 5 Hz). EDA Guidelines: generally < 10 Hz |
| **2. Artifact removal** | Threshold criterion | Reject samples where change exceeds threshold (e.g., > 0.01 μS/sample) |
| **3. (Optional) Decomposition** | Tonic vs phasic | cvxEDA, highpass (0.05 Hz), or median smoothing |

### NeuroKit2 usage
```python
import neurokit2 as nk
eda_clean = nk.eda_clean(eda_signal, sampling_rate=100, method="neurokit")  # 3 Hz lowpass
# Or full pipeline:
signals, info = nk.eda_process(eda_signal, sampling_rate=100)
```

### Automated quality (large datasets)
- **[EDAQA](https://github.com/iankleckner/EDAQA)**: Detects out-of-range, fast changes, wear status, invalid segments

---

## RSA (Respiratory Sinus Arrhythmia)

### Noise sources
- **Technical**: Power line, baseline wander, muscle artifacts
- **Physiological**: Ectopic beats (abnormal RR intervals)

### Standard preprocessing

**If RSA is a continuous waveform** (your case — already derived from HR/HRV):
- **Bandpass filter**: Keep **0.12–0.4 Hz** (HF band, respiratory range)
- Removes DC drift, high-freq noise

**If starting from RR intervals**:
- Remove ectopic beats (e.g., Kubios, NeuroKit2 `correct_artifacts=True`)
- Interpolate to regular grid before spectral analysis
- Two-step: CEEMDAN for Gaussian noise, adaptive filter for ectopics ([Nature 2022](https://www.nature.com/articles/s41598-022-21776-2))

### NeuroKit2 HRV bands (for reference)
- HF (RSA): **0.15–0.4 Hz**
- LF: 0.04–0.15 Hz
- VLF: 0.0033–0.04 Hz

### For continuous RSA waveform at 100 Hz
```python
# Bandpass 0.12–0.4 Hz to isolate respiratory modulation
rsa_clean = nk.signal_filter(rsa_signal, sampling_rate=100, lowcut=0.12, highcut=0.4, order=4)
```

---

## Open Source Tools

| Tool | Signals | Key functions |
|------|---------|---------------|
| **[NeuroKit2](https://neuropsychology.github.io/NeuroKit/)** | EDA, ECG, PPG, RSP, HRV | `eda_clean`, `eda_process`, `signal_filter`, `ecg_clean`, `hrv()` |
| **[pyEDA](https://github.com/akaraspt/pyEDA)** | EDA | Preprocessing + feature extraction |
| **[PhysioKit](https://ambiqai.github.io/physiokit/)** | ECG, PPG, RSP, HRV | Wearable-focused, shared filtering utils |
| **[EDAQA](https://github.com/iankleckner/EDAQA)** | EDA | Automated quality assessment |

---

## Recommended pipeline for your EDA + RSA data (100 Hz)

1. **EDA**: `nk.eda_clean(eda, sampling_rate=100, method="neurokit")` — 3 Hz lowpass
2. **RSA**: `nk.signal_filter(rsa, sampling_rate=100, lowcut=0.12, highcut=0.4)` — bandpass
3. (Optional) Remove motion artifacts via threshold or EDAQA
4. Run your overlap/covariance analysis on cleaned signals

---

## References

- EDA Guidelines: https://edaguidelines.github.io/analysis/pre
- NeuroKit2 EDA: https://neuropsychology.github.io/NeuroKit/functions/eda.html
- NeuroKit2 signal_filter: https://neuropsychology.github.io/NeuroKit/functions/signal.html
- HRV HF band (RSA): 0.12–0.4 Hz resting (MindWare, Task Force 1996)
- Nature 2022: Two-step HRV denoising (Gaussian + ectopic)
