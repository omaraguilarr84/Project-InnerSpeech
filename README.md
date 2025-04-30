# 🧠 Lexical EEG Decoding: Spatial Dynamics & Functional Connectivity

**Paper Title**: *Spatial Dynamics and Functional Connectivity in EEG: Insights from Lexical Processing*  
**Authors**: Omar Aguilar, Ramon Grullon  
**Affiliation**: Georgia Institute of Technology, Atlanta, GA, USA  
**DOI**: [10.1101/2025.04.29.651245](https://doi.org/10.1101/2025.04.29.651245)

---

## 📚 Overview

This repository accompanies our research on lexical processing using EEG data. Our study investigates the **spatial dynamics** and **functional connectivity** of the brain during silent reading tasks using **Random Forest classifiers**. We demonstrate how combining traditional frequency-domain EEG features with machine learning enables high-accuracy classification of semantic word categories (social vs. numeric).

---

## 🔬 Abstract

Lexical processing requires semantic integration across distributed neural networks. Using EEG recordings from a silent reading task, we extracted band power and coherence features to classify word categories. Our **spatial dynamics model** achieved **100% accuracy**, identifying theta and alpha band activity in parietal and occipital regions as key. A complementary **functional connectivity model** reached **85% accuracy**, revealing the significance of inter-regional coherence. These results lay foundational work for future EEG-based brain-computer interfaces (BCIs).

---

## 🧪 Methods

### 🧠 EEG Dataset
- Source: [Liwicki et al., 2023 Scientific Data](https://www.nature.com/articles/s41597-023-02249-4)
- Participants: 5 total (1 used for analysis due to data quality)
- Trials: 320 trials per participant
- Task: Silent inner-speech repetition of 8 words (4 social, 4 numeric)

### ⚙️ Feature Extraction
- **Frequency Bands**:  
  - Theta: 4–8 Hz  
  - Alpha: 8–13 Hz  
  - Beta: 13–30 Hz  
  - Gamma: 30–50 Hz
- **Spatial Features**: Band power from 25 electrodes (parietal, occipital, frontal, etc.)
- **Connectivity Features**: Coherence between filtered electrode pairs

### 🤖 Models
- **Model 1**: Band Power → Random Forest → 100% accuracy
- **Model 2**: Coherence → Random Forest → 85% accuracy
- Grid search with 5-fold cross-validation for hyperparameter tuning

---

## 📊 Results

### ✅ Spatial Dynamics Model
- Perfect classification of word categories (social vs. numeric)
- Theta and alpha bands from parietal and occipital electrodes showed highest importance

### 🔗 Functional Connectivity Model
- Revealed strong coherence between parietal and central electrodes
- Unique isolated connections (e.g., AF8–PO7) suggest specialized processing pathways

---

## 💡 Implications

- Demonstrates feasibility of EEG-based semantic classification
- Enables design of compact BCI systems using minimal electrode coverage
- Bridges traditional EEG analysis with interpretable machine learning
