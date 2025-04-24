#  Order of Preprocessing and Analysis Code

1. **`HBN_preprocessing_resting`**  
   - Loads EEG files  
   - Applies notch filtering, bandpass filtering, and bad channel removal  
   - 📌 *Refer to the Healthy Brain Network guidelines for the data download process.*

2. **`HBN_resting_event_extract_rp`**  
   - Extracts segments based on event logs:  
     - Eyes-open: 16.8 seconds  
     - Eyes-closed: 33.6 seconds  
   - Computes relative phase from the extracted segments

3. **`HBN_kmean_4mean_combined_centroid`**  
   - Performs 4-means clustering on relative phase data  
   - Uses a predefined **combined centroid** as the clustering reference  
   - 🔍 *See the `library/` folder for centroid and related functions*

4. **`HBN_calculate_beta_by_subject`**  
   - Computes **beta values** as described in the paper  
   - Calculates additional subject-wise metrics such as:
     - Mode ratio
     - Dwell time

5. **`HBN_tbr_cal`**  
   - Calculates the **Theta/Beta Ratio (TBR)** from the EEG data

6. **`HBN_logistic_roc`**  
   - Combines relative phase-based variables using a **logistic function**  
   - Plots **ROC curves** to evaluate classification performance

7. **`HBN_logistic_corr`**  
   - Combines variables using a logistic function  
   - Computes correlation with the **SWAN inattentive score**  
     *(SWAN = Strengths and Weaknesses of ADHD Symptoms and Normal Behavior Scale)*
