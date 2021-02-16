# Translation Initiation Site Prediction in _Arabidopsis thaliana_ Using Synthetic Datasets and Black-Box Models

### Bachelor's dissertation in Center of Biotech Data Science, Ghent University Global Campus
#### Yunseol Park
#### Counselors: [Espoir Kabanga](https://github.com/EspoirKabanga), Jasper Zuallaert
#### Supervisors: [Arnout Van Messem](https://github.com/avmessem), Wesley De Neve

## Description
The main goals of the research effort presented in this project are as follows:
- Identify meaningful features of the TIS prediction model.
- Compare the true black-box model with the synthetic black-box model.
- Investigate the effect of noisy data on the TIS prediction.

In order to achieve these goals, the following steps are taken:
1. Generate the TIS synthetic dataset.
2. Train the model with real and synthetic data (_A. thaliana_).
3. Compare the results of the models trained on real and synthetic data.
4. Perform feature analysis.
5. Train the prediction model with noisy data.

<p align="center">
  <img src="https://github.com/YunseolPark/BAthesis-TIS-prediction-using-synthetic-data/blob/master/images/schematic_plan.png" />
  Schematic Plan
</p>

*The following readme file is not a complete account of the bachelor's dissertation. It is merely the explanation of the python files uploaded to this page with regards to the dissertation.*

## Table of Contents
1. [Introduction](#introduction)
2. [Dataset](#dataset)
3. [Model Training](#model-training)
4. [Feature Analysis](#feature-analysis)
5. [Noise Analysis](#noise-analysis)

## Introduction

Prediction of translation initiation sites (TISs) can give insight into translation and the proteins synthesized by certain mRNAs. Thus, it is important for genome analysis and annotation. Furthermore, the mechanisms of translation have not been perfectly studied. Therefore, by interpreting the prediction model, it may even aid in uncovering new translation mechanisms or give emphasis to an existing one.
However, a lot of real-world datasets contain noise and many genome annotations like TIS prediction are high-risk problems. Furthermore, real-world data are extremely complex, which makes it difficult to find the features that influence the decision of the model.
Synthetic data can be used to solve this problem. In particular, synthetic data can be used to give insight into the features of the model and thus into the real-world data. They are suitable for this purpose since they are constructed by selecting and incorporating some of the complex features of the real-world dataset. The outcome of the synthetic model can then be compared to the real model to find the features that contribute most to the prediction of TIS. Furthermore, the effect of noise on datasets can also be investigated to see how the model performs with noisy data.

## Dataset

### Real Dataset
The real dataset used was a dataset of Arabidopsis thaliana TIS dataset taken from [Magana-Mora et al. (2013)]. It is a balanced dataset, with both positive and negative sets containing 207102 sequences that are 300 nucleotides long.
The dataset is constructed in a way that the TIS (positive set) or non-TIS ATG (negative set) are centrally located, with 148 sequences upstream and 149 sequences downstream.

### Synthetic Dataset
The synthetic dataset, generated via the Python code in file `GenerateTIS.py`, contains 27102 sequences that are 300 nucleotides long. Each TIS is centrally located in the sequence, on the 150th – 152nd nucleotide. They were generated with the same structure as the real dataset.
The dataset is generated with 5 features, adding each feature in a different step (Figure 1). The following sections will discuss each step in detail.

<p align="center">
  <img src="https://github.com/YunseolPark/BAthesis-TIS-prediction-using-synthetic-data/blob/master/images/synthetic_data.png" />
  Process used for generating a synthetic dataset, containing positive (left) and negative (right) samples.
  The central start codon is in bold and underlined, while the changes are indicated in red.
</p>

#### Basic Structure
The basic structure is designed so that each TIS, ‘ATG’, is centrally located, with 150 nucleotides upstream and 150 nucleotides downstream, indicated as ‘u’s and ‘d’s in Figure 6 (1), respectively. The upstream and downstream sequences are assigned equal lengths which are divisible by 3. The sequence is designed in this way for easy manipulation, namely, easier addition of codons. Then, after the dataset has been created, sequences on the sides will be cleaved off to get the desired length of 300.

#### Consensus (Kozak) Sequence
The first feature to be added to the synthetic dataset is the consensus sequence, as shown in Figure 6 (2). The upstream consensus sequence spans from positions -10 to -1 and the downstream consensus sequence spans from positions +3 to +12 where 0, 1, 2 indicate the ‘ATG’ (Table 3). In this dissertation, a nucleotide followed by a superscript will be used to denote a nucleotide at that position. For example, G at position +3 will be denoted as G_+3.
The data for the consensus sequence was obtained from the real dataset using the Python code,`ConsensusSequence.py`. The consensus sequence was determined using the 50/70% rule, a modification of Cavener 50/75% rule ([Cavener, 1987]).

|           | -10 | -9 | -8 | -7 | -6 | -5 | -4 |  -3 |  -2 | -1 |
|:---------:|:---:|:--:|:--:|:--:|:--:|:--:|:--:|:---:|:---:|:--:|
|     A%    |  35 | 32 | 34 | 35 | 36 | 33 | 45 |  50 |  42 | 44 |
|     C%    |  16 | 17 | 20 | 17 | 16 | 24 | 14 |  11 |  29 | 19 |
|     G%    |  20 | 21 | 18 | 20 | 22 | 17 | 21 |  24 |  9  | 23 |
|     T%    |  29 | 30 | 27 | 28 | 27 | 26 | 19 |  15 |  20 | 14 |
| Consensus |  a  |  a |  a |  a |  a |  a |  a | A/G | A/C |  a |

<br/>

|           | -10 | -9 | -8 | -7 | -6 | -5 | -4 |  -3 |  -2 | -1 |
|:---------:|:---:|:--:|:--:|:--:|:--:|:--:|:--:|:---:|:---:|:--:|
|     A%    |  35 | 32 | 34 | 35 | 36 | 33 | 45 |  50 |  42 | 44 |
|     C%    |  16 | 17 | 20 | 17 | 16 | 24 | 14 |  11 |  29 | 19 |
|     G%    |  20 | 21 | 18 | 20 | 22 | 17 | 21 |  24 |  9  | 23 |
|     T%    |  29 | 30 | 27 | 28 | 27 | 26 | 19 |  15 |  20 | 14 |
| Consensus |  a  |  a |  a |  a |  a |  a |  a | A/G | A/C |  a |


#### Upstream ATG
An upstream ATG is only inserted for the positive dataset, as shown in Figure 6 (3). It is added to any random position between the start of the sequence and the start of the consensus sequence. The number of start codons to be inserted is determined randomly, varying from 0 to 2.
[Zuallaert et al. (2018b)] reported that upstream ATGs have a negative influence on the prediction of TIS while downstream ATGs have a positive influence by learning the properties of the scanning model. However, in this dissertation, it was added as a feature for the positive set in order to investigate if the model would be able to learn the leaky scanning hypothesis along with the scanning model given enough context around TIS.
Furthermore, it was added to investigate if the positive influence of the context around TIS would be able to compensate for the negative influence of upstream ATGs. Since almost 40% of TISs from the real-world data have upstream ATGs ([Pedersen & Nielsen, 1997]), this investigation would give insight into how the prediction model behaves with TISs having upstream ATGs.

#### Downstream Stop Codons
The downstream stop codons are only inserted for the negative dataset as it signals the end of the translation mechanism. This can be seen in Figure 6 (3). One of the stop codons, ‘TAA’, ‘TAG’, or ‘TGA’, is randomly selected and added at a random position between the end of the central ATG and the end of the sequence.
As this feature is only added to the negative dataset, any in-frame downstream stop codons added to the positive dataset due to the random addition of nucleotides are removed from the sequence.

#### Donor Splice Site
The donor splice site pattern, as seen in Table 4, shows the first three consensus nucleotide frequencies of donor splice site. These three nucleotides are a part of the splicing residues that can be found in the exon.
The donor splice site is added downstream for the samples in the positive dataset and upstream for the samples in the negative dataset, as illustrated in Figure 6 (4). This is because, normally, the first gene (exon) will contain a TIS while non-start ATGs are more likely to occur inside a gene ([Zuallaert et al., 2018b]).

|  Position |  -3 |  -2 |  -1 |
|:---------:|:---:|:---:|:---:|
|     A     | 352 | 651 |  85 |
|     C     | 361 | 131 |  45 |
|     G     | 154 |  87 | 777 |
|     T     | 132 | 159 |  91 |
| Consensus |  C  |  A  |  G  |

#### Nucleotide Frequency
The empty space after the insertion of codons is then filled with nucleotide frequency, as seen in Figure 6 (5). The nucleotide frequency for the samples in the positive dataset and negative dataset can be found below.  can be found in Table 5 and the nucleotide frequency for the samples in the negative dataset can be found in Table 6.

| Class |  A%  |  C%  |  G%  |  T%  |
|:-----:|:----:|:----:|:----:|:----:|
|  UTR  | 31.1 | 19.6 | 15.5 | 33.8 |
|  ORF  | 26.1 | 23.0 | 20.4 | 29.7 |

<br/>

|    Class   |  A%  |  C%  |  G%  |  T%  |
|:----------:|:----:|:----:|:----:|:----:|
|  Upstream  | 31.4 | 18.4 | 18.5 | 31.7 |
| Downstream | 31.3 | 18.1 | 18.5 | 31.4 |

## Model Training and Evaluation
Models are trained using both synthetic and real datasets via TISRover API ([Zuallaert et al., 2018b]). The synthetic black-box model (SBBM) is a model trained on synthetic data only while the real black-box model (RBBM) is a model trained on real data only. The combined black-box model (CBBM) is a model trained on a combination of synthetic and real data in a 1:1 ratio.

All three models were trained with datasets partitioned in a way that 90% is for training and validation sets (further partitioned into 7/8th and 1/8th respectively), and 10% is for test set.

### Metrics
Recall (sensitivity), precision, accuracy, and F1 score were used for the evaluation metrics.

![Recall](https://github.com/YunseolPark/BAthesis-TIS-prediction-using-synthetic-data/blob/master/images/formula_recall.png)
<br/>
![precision](https://github.com/YunseolPark/BAthesis-TIS-prediction-using-synthetic-data/blob/master/images/formula_precision.png)
<br/>
![accuracy](https://github.com/YunseolPark/BAthesis-TIS-prediction-using-synthetic-data/blob/master/images/formula_accuracy.png)
<br/>
![F1](https://github.com/YunseolPark/BAthesis-TIS-prediction-using-synthetic-data/blob/master/images/formula_f1.png)

where TP, TN, FP, FN refer to true positive, true negative, false positive, and false negative, respectively.

## Feature Analysis

### Perturbance-based Method
Perturbation (Shrikumar et al., 2017) is the modification of an input feature where the modification can be used to analyze the impact of said feature. The perturbance strategy used in this project is occlusion, which is the removal of feature(s)/region(s) in the input. It is possible to occlude a single feature and it is also possible to occlude all features except for the feature of interest.

<p align="center">
  <img src="https://github.com/YunseolPark/BAthesis-TIS-prediction-using-synthetic-data/blob/master/images/occlusion.png" />
</p>

The analysis of models with missing features allows for the evaluation and visualization of important features through comparison of the evaluation metrics (obtained by testing with the real dataset). The lower the metrics value, the more influential the (missing) feature and vice versa.

The features and models to be analyzed are:
- consensus sequence: SBBM<sub>Consensus</sub>
- upstream ATG: SBBM<sub>ATG</sub>
- downstream stop codon: SBBM<sub>Stop</sub>
- donor splice site: SBBM<sub>Donor</sub>

The analysis of models with single features also allows for the evaluation and visualization of important features through training SBBM with only one feature at a time. The higher the evaluation metrics (obtained by testing with the real dataset), the more influential the (single) feature and vice versa.

The features to be analyzed are:
- nucleotide frequency
- consensus sequence
- upstream ATG
- downstream stop codon
- donor splice site
- nucleotide triplet

The features consensus sequence, upstream ATG, downstream stop codon, and donor splice site also contain the feature nucleotide frequency as the base in generating the synthetic dataset. The nucleotide triplet is added to test the impact of codons on TIS prediction.

<table>


### Attribution-based Method

#### Sequence Logo
The sequence logo can visualize the contribution of each nucleotide in a certain position. The contribution scores are processed so that on a certain position, averages of each nucleotide are obtained for true positive samples. Then it is visualized on a sequence logo created from Logomaker. The code for obtaining the averages can be found in `ContributionScore.py` and `LogoFile.py`, and the plotting of the averages can be found in `PlotLogo.py`.

## Noise Analysis

Noisy data can affect a model by reducing its performance in general. Noise can be measured by equalized loss of accuracy (ELA). It measures the behavior of a model at a certain noise level (Sáez et al., 2016).

![ELA](https://github.com/YunseolPark/BAthesis-TIS-prediction-using-synthetic-data/blob/master/images/formula_ela.png)

where [𝐸𝐿𝐴]<sub>𝑥%</sub> refers to the measure of behavior at noise level 𝑥%, 𝐴<sub>0%</sub> refers to the accuracy without any noise, and 𝐴<sub>𝑥%</sub> refers to the accuracy at noise level 𝑥%. Accuracy is calculated  and then converted into a percentage value. A lower ELA value will indicate that the model behaves better with the given noise.

There are many sources of noise, but in this dissertation, only two are investigated:
1. incorrect labeling or misclassification of training data (class noise); and
2. incorrect features (attribute noise).

Class noise, more specifically misclassification, is investigated via introducing labeling errors. The effects of noise during training a model is investigated by introducing 10% (𝑥= 10) noise and 20% (𝑥= 20) noise artificially into the datasets.

Attribute noise is investigated via introducing ‘incorrect’ features. In this case, these features are not incorrect, but are simply features that were added to the synthetic dataset. Attribute noise is difficult to investigate using real datasets and is almost impossible to induce artificially. Thus, only the synthetic dataset is used.
Each model that was trained with a dataset having a missing feature is used as a ‘noiseless’ model and ‘noise’ is added by adding the dataset with all features (Figure 8). They are then tested on the real test set. Through this, it is anticipated that the results from Section 3.5 will be validated. If the feature is important, it will not be considered as noise, and the ELA value will decrease. Conversely, if the feature is not important, then the ELA value will increase.

During this experiment, the codon usage table is also considered. The codon usage table is a table that shows the frequency of each codon. The code for obtaining codon usage can be found in `CodonUsage.py`. The feature of codons will be denoted as a nucleotide triplet in this documentation.
By investigating the impact of ‘noise’ that may be added with the insertion of nucleotide triplets, the impact of codons on TIS prediction can be investigated. The dataset with nucleotide triplets is created by removing nucleotide frequency and replacing it with the values obtained from the codon usage table. This would eliminate the addition of nonsense codons and will ensure that the dataset has meaningful triplets, although it would not ensure the production of functional proteins from the sequences.


<References>

[Magana-Mora et al. (2013)]: <https://academic.oup.com/bioinformatics/article/29/1/117/272605>
[Cavener, 1987]: <https://academic.oup.com/nar/article-abstract/15/4/1353/2377965?redirectedFrom=fulltext>
[Zuallaert et al. (2018b)]: <https://dl.acm.org/doi/10.5555/3282702.3282707>
[Zuallaert et al., 2018b]: <https://dl.acm.org/doi/10.5555/3282702.3282707>
[Pedersen & Nielsen, 1997]: <https://pubmed.ncbi.nlm.nih.gov/9322041/>
