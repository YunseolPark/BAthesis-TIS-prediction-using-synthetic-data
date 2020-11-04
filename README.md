# Translation Initiation Site Prediction in Arabidopsis thaliana Using Synthetic Datasets and Black-Box Models

#### Yunseol Park
#### Supervisors: [Espoir Kabanga](https://github.com/EspoirKabanga), Jasper Zuallaert, [Arnout Van Messem](https://github.com/avmessem), Wesley De Neve

## Description
The main goals of the project are to create a translation initiation site prediction model and investigate the black-box nature of the model using synthetic datasets.
The synthetic datasets were built from scratch using features known to be important in both biological and machine learning aspects.

## Table of Contents
1. [Introduction](#introduction)
2. [Dataset](#dataset)
3. [Model Training](#model-training)
4. [Feature Analysis](#feature-analysis)
5. [Noise Analysis](#noise-analysis)
6. [References](#references)

## Introduction

Prediction of translation initiation sites (TISs) can give insight into translation and the proteins synthesized by certain mRNAs. Thus, it is important for genome analysis and annotation. Furthermore, the mechanisms of translation have not been perfectly studied. Therefore, by interpreting the prediction model, it may even aid in uncovering new translation mechanisms or give emphasis to an existing one.
However, a lot of real-world datasets contain noise and errors and many genome annotations like TIS prediction are high-risk problems. Furthermore, real-world data are extremely complex, which makes it difficult to find the features that influence the decision of the model.
Synthetic data can be used to solve this problem. In particular, synthetic data can be used to give insight into the features of the model and thus into the real-world data. They are suitable for this purpose since they are constructed by selecting and incorporating some of the complex features of the real-world dataset. The outcome of the synthetic model can then be compared to the real model to find the features that contribute most to the prediction of TIS. Furthermore, the effect of noise on datasets can also be investigated to see how the model performs with noisy data.

The main goals of the research effort presented in this project are as follows:
- Identify meaningful features of the TIS prediction model.
- Compare the true black-box model with the synthetic black-box model.
- Investigate the effect of noisy data on the TIS prediction.

In order to achieve these goals, the following steps are taken:
1. Generate the TIS synthetic dataset.
2. Train the model with real and synthetic data (A. thaliana).
3. Compare the results of the models trained on real and synthetic data.
4. Perform feature analysis.
5. Train the prediction model with noisy data.

## Dataset

### Real Dataset
The real dataset used was a dataset of Arabidopsis thaliana TIS dataset taken from Magana-Mora et al. (2013). It is a balanced dataset, with both positive and negative sets containing 207102 sequences that are 300 nucleotides long.
The dataset is constructed in a way that the TIS (positive set) or non-TIS ATG (negative set) are centrally located, with 148 sequences upstream and 149 sequences downstream.

### Synthetic Dataset
The synthetic dataset, generated via the Python code in file `GenerateTIS.py`, contains 27102 sequences that are 300 nucleotides long. Each TIS is centrally located in the sequence, on the 150th – 152nd nucleotide. They were generated with the same structure as the real dataset.
The dataset is generated with 5 features, adding each feature in a different step (Figure 6). The following sections will discuss each step in detail.

<image>

#### Basic Structure

The basic structure is designed so that each TIS, ‘ATG’, is centrally located, with 150 nucleotides upstream and 150 nucleotides downstream, indicated as ‘u’s and ‘d’s in Figure 6 (1), respectively. The upstream and downstream sequences are assigned equal lengths which are divisible by 3. The sequence is designed in this way for easy manipulation. The sequences’ divisibility allows easier addition of codons.
In designing the sequence, only the canonical TIS was considered for this dissertation.

#### Consensus (Kozak) Sequence

The first feature to be added to the synthetic dataset is the consensus sequence, as shown in Figure 6 (2). It is only added to the positive dataset since it is a sequence that would indicate the presence of TIS. The upstream consensus sequence spans from positions -10 to -1 and the downstream consensus sequence spans from positions +3 to +12 where 0, 1, 2 indicate the ‘ATG’ (Table 3). In this dissertation, a nucleotide followed by a superscript will be used to denote a nucleotide at that position. For example, G at position +3 will be denoted as G+3.
The data for the consensus sequence was obtained from the real dataset using the Python code, `ConsensusSequence.py`. The consensus sequence was determined using the 50/70% rule, a modification of Cavener 50/75% rule.

|           | -10 | -9 | -8 | -7 | -6 | -5 | -4 |  -3 |  -2 | -1 |
|:---------:|:---:|:--:|:--:|:--:|:--:|:--:|:--:|:---:|:---:|:--:|
|     A%    |  35 | 32 | 34 | 35 | 36 | 33 | 45 |  50 |  42 | 44 |
|     C%    |  16 | 17 | 20 | 17 | 16 | 24 | 14 |  11 |  29 | 19 |
|     G%    |  20 | 21 | 18 | 20 | 22 | 17 | 21 |  24 |  9  | 23 |
|     T%    |  29 | 30 | 27 | 28 | 27 | 26 | 19 |  15 |  20 | 14 |
| Consensus |  a  |  a |  a |  a |  a |  a |  a | A/G | A/C |  a |


|           | -10 | -9 | -8 | -7 | -6 | -5 | -4 |  -3 |  -2 | -1 |
|:---------:|:---:|:--:|:--:|:--:|:--:|:--:|:--:|:---:|:---:|:--:|
|     A%    |  35 | 32 | 34 | 35 | 36 | 33 | 45 |  50 |  42 | 44 |
|     C%    |  16 | 17 | 20 | 17 | 16 | 24 | 14 |  11 |  29 | 19 |
|     G%    |  20 | 21 | 18 | 20 | 22 | 17 | 21 |  24 |  9  | 23 |
|     T%    |  29 | 30 | 27 | 28 | 27 | 26 | 19 |  15 |  20 | 14 |
| Consensus |  a  |  a |  a |  a |  a |  a |  a | A/G | A/C |  a |



#### Upstream ATG

An upstream ATG is only inserted for the positive dataset, as shown in Figure 6 (3). It is added to any random position between the start of the sequence and the start of the consensus sequence. The number of start codons to be inserted is determined randomly, varying from 0 to 2.
Zuallaert et al. (2018b) reported that upstream ATGs have a negative influence on the prediction of TIS while downstream ATGs have a positive influence by learning the properties of the scanning model. However, in this dissertation, it was added as a feature for the positive set in order to investigate if the model would be able to learn the leaky scanning hypothesis along with the scanning model given enough context around TIS.
Furthermore, it was added to investigate if the positive influence of the context around TIS would be able to compensate for the negative influence of upstream ATGs. Since almost 40% of TISs from the real-world data have upstream ATGs (Pedersen & Nielsen, 1997), this investigation would give insight into how the prediction model behaves with TISs having upstream ATGs.

#### Downstream Stop Codons

The downstream stop codons are only inserted for the negative dataset as it signals the end of the translation mechanism. This can be seen in Figure 6 (3). The stop codons are ‘TAA’, ‘TAG’, and ‘TGA’. One stop codon is randomly selected and added at a random position between the end of the central ATG and the end of the sequence.
As this feature is only added to the negative dataset, any in-frame downstream stop codons added to the positive dataset due to the random addition of nucleotides are removed from the sequence.

#### Donor Splice Site

The donor splice site pattern, as seen in Table 4, shows the first three consensus nucleotide frequencies of donor splice site. These three nucleotides are a part of the splicing residues that can be found in the exon.
The donor splice site is added downstream for the samples in the positive dataset and upstream for the samples in the negative dataset, as illustrated in Figure 6 (4). This is because, normally, the first gene (exon) will contain a TIS while non-start ATGs are more likely to occur inside a gene (Zuallaert, Kim, et al., 2018).
The data for the donor splice site consensus sequence were obtained from Kim (2019).

|  Position |  -3 |  -2 |  -1 |
|:---------:|:---:|:---:|:---:|
|     A     | 352 | 651 |  85 |
|     C     | 361 | 131 |  45 |
|     G     | 154 |  87 | 777 |
|     T     | 132 | 159 |  91 |
| Consensus |  C  |  A  |  G  |


#### Nucleotide Frequency

The empty space after the insertion of codons is then filled with nucleotide frequency, as seen in Figure 6 (5). The nucleotide frequency for the samples in the positive dataset can be found in Table 5 and the nucleotide frequency for the samples in the negative dataset can be found in Table 6.


| Class |  A%  |  C%  |  G%  |  T%  |
|:-----:|:----:|:----:|:----:|:----:|
|  UTR  | 31.1 | 19.6 | 15.5 | 33.8 |
|  ORF  | 26.1 | 23.0 | 20.4 | 29.7 |

|    Class   |  A%  |  C%  |  G%  |  T%  |
|:----------:|:----:|:----:|:----:|:----:|
|  Upstream  | 31.4 | 18.4 | 18.5 | 31.7 |
| Downstream | 31.3 | 18.1 | 18.5 | 31.4 |
 
The nucleotide frequency was gathered from the real dataset. Only the regions disregarding the central ATG codon and the consensus sequence were taken for nucleotide frequency. The reason behind this was to use nucleotide frequency to replace the consensus sequence when training a model without it in a later part of the experiment.

## Model Training

Models are trained using both synthetic and real datasets. The synthetic black-box model (SBBM) is a model trained on synthetic data only while the real black-box model (RBBM) is a model trained on real data only. The results obtained by the two models will be compared.

The combined black-box model (CBBM) is a model trained on a combination of synthetic and real data. The length of the input data is set to be the same for both the synthetic dataset and the real dataset, and they are added in a 1:1 ratio.
CBBM can be used to compare the results of SBBM and RBBM. Since the two will most likely learn differently, combining them may allow the model to learn features from both datasets. If this is the case, then CBBM will have similar results to that of RBBM. However, if SBBM and RBBM learn from different features, then CBBM will be different from RBBM. In this case, RBBM may have better results than the two by learning both features, or, as it was for the case of the experiment by Kim (2019), it may have worse results if SBBM learns from insufficient features.

Recall (Sensitivity), Precision, Accuracy, and F1 score were used for the evaluation metrics.

## Feature Analysis

### Perturbance-based Method
The impact of features in TIS prediction can be investigated by training a model with synthetic datasets having modified features. This is a perturbance-based strategy (Shrikumar et al., 2017). Perturbation is the modification of an input feature where the modification can be used to analyze the impact of said feature. Many strategies could be implemented using perturbation, but in this dissertation, occlusion is used. Occlusion is the removal of feature(s) or region(s) in the input. It is possible to occlude a single feature and it is also possible to occlude all features except for the feature of interest. From the results, the impact of the feature can be analyzed on the prediction model. In this dissertation, two types of models with occlusion are analyzed: models with missing features, and models with single features.
The analysis of models with missing features allows for the evaluation and visualization of important features through comparison of the evaluation metrics. The metrics are obtained by testing with the real dataset.
The models that are trained on these data will be denoted by SBBM, followed by a subscript of the missing feature. The features to be analyzed are (1) consensus sequence (SBBMConsensus), (2) upstream ATG (SBBMATG), (3) downstream stop codon (SBBMStop), and (4) donor splice site (SBBMDonor).
The analysis will find the features that have an impact on prediction. If the evaluation metrics without a feature gives a lower value than that of normal SBBM, then it can be concluded that the feature has a positive influence on the prediction of TIS. The lower the value, the more influential the (missing) feature and vice versa.
The analysis of models with missing features allows for the evaluation and visualization of important features through comparison of the evaluation metrics. The metrics are obtained by testing with the real dataset.

The models that are trained on these data will be denoted by SBBM, followed by a subscript of the missing feature. The features to be analyzed are (1) consensus sequence (SBBMConsensus), (2) upstream ATG (SBBMATG), (3) downstream stop codon (SBBMStop), and (4) donor splice site (SBBMDonor).
The analysis will find the features that have an impact on prediction. If the evaluation metrics without a feature gives a lower value than that of normal SBBM, then it can be concluded that the feature has a positive influence on the prediction of TIS. The lower the value, the more influential the (missing) feature and vice versa.
The analysis of models with single features also allows for the evaluation and visualization of important features. The evaluation of the features is to be performed by training SBBM with only one feature at a time and testing it on the real test set. The features to be analyzed are (1) nucleotide frequency, (2) consensus sequence, (3) upstream ATG, (4) downstream stop codon, (5) donor splice site, and (6) nucleotide triplet. The features consensus sequence, upstream ATG, downstream stop codon, and donor splice site also contain the feature nucleotide frequency as the base in generating the synthetic dataset. The nucleotide triplet is added to test the impact of codons on TIS prediction. The details on this feature are further explained in Section 3.6.
The analysis will further validate the features that impact the prediction of TIS. Since the nucleotide frequency is used as the base for all other features, if the evaluation metrics of a model with a single feature has a higher value than that of nucleotide frequency, then it can be concluded that it has a positive influence on prediction. The higher the value, the more influential the (single) feature and vice versa.

### Atribution-based Method

#### Sequence Logo
The sequence logo can visualize the contribution of each nucleotide in a certain position. The contribution scores are processed so that on a certain position, averages of each nucleotide are obtained for true positive samples. Then it is visualized on a sequence logo created from Logomaker. The code for obtaining the averages can be found in Appendix C.

## Noise Analysis

Noisy data can affect a model by reducing its performance in general. There are many sources of noise, but in this dissertation, only two are investigated:
1. incorrect labeling or misclassification of training data (class noise); and
2. incorrect features (attribute noise).
Noise can be measured by equalized loss of accuracy (ELA). It measures the behavior of a model at a certain noise level (Sáez et al., 2016).
𝐸𝐿𝐴𝑥%= 100−𝐴𝑥%𝐴0%
(5)
where 𝐸𝐿𝐴𝑥% refers to the measure of behavior at noise level 𝑥%, 𝐴0% refers to the accuracy without any noise, and 𝐴𝑥% refers to the accuracy at noise level 𝑥%. Accuracy is calculated by Equation (3) and then converted into a percentage value. In this dissertation, an ELA between noise levels 𝑎 and 𝑏 will be indicated as 𝐸𝐿𝐴𝑎−𝑏%, where 𝑎<𝑏. In the interpretation of ELA values, a lower value will indicate that the model behaves better with the given noise.
ELA is a measurement that considers both the noise robustness and the prediction rate of a model with the introduction of certain noise (Sáez et al., 2016), where robustness refers to the stability of a prediction model when noise is introduced. This characteristic is important since many other measurements consider only one of two, which makes them not very conclusive as the two concepts may imply two different conclusions. Furthermore, it takes into account a model’s behavior with clean data. Thus, ELA is suitable for comparing behavior with noise in different models (Nazari et al., 2018).

Class noise, more specifically misclassification, is investigated via introducing labeling errors. The effects of noise during training a model is investigated by introducing 10% (𝑥= 10) noise and 20% (𝑥= 20) noise artificially into the datasets. During the partitioning of data, 10% of the training set will be switched between the positive and negative datasets while the labels will be left as they are to imitate labeling errors (Figure 7). This is done to both the real and synthetic datasets.

Attribute noise is investigated via introducing ‘incorrect’ features. In this case, these features are not incorrect, but are simply features that were added to the synthetic dataset. Attribute noise is difficult to investigate using real datasets and is almost impossible to induce artificially. Thus, only the synthetic dataset is used.

Each model that was trained with a dataset having a missing feature (from Section 3.5) is used as a ‘noiseless’ model and ‘noise’ is added by adding the dataset with all features (Figure 8). They are then tested on the real test set. Through this, it is anticipated that the results from Section 3.5 will be validated. If the feature is important, it will not be considered as noise, and the ELA value will decrease with the addition of the feature in the dataset. Conversely, if the feature is not important, then the ELA value will increase.
During this experiment, the codon usage table is also considered, aside from the datasets used in Section 3.5. The codon usage table is a table that shows the frequency of each codon and can be found in Appendix D. The feature of codons will be denoted as a nucleotide triplet for the rest of this dissertation.
A nucleotide triplet is considered as a feature to check if the use of codons would have an impact on TIS prediction. It was mentioned in many articles that these codons are factors that have an influence on genome prediction. Kim (2019) also speculated that disregard for the relationship of nucleotides, and thus the production of nonsense codons, may have been one of the reasons for the poor prediction rate in SBBM.
By investigating the impact of ‘noise’ that may be added with the insertion of nucleotide triplets, the impact of codons on TIS prediction can be investigated. The dataset with nucleotide triplets is created by removing nucleotide frequency and replacing it with the values obtained from the codon usage table. This would eliminate the addition of nonsense codons and will ensure that the dataset has meaningful triplets, although it would not ensure the production of functional proteins from the sequences.

## References
