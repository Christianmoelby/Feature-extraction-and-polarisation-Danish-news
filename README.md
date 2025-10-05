# Feature-extraction-Danish-news

Contains essential code to scrape and collect Danish news and extract features from text, images, and faces. Formats features into a token-based dataset.

Extended with practical implementation of Maximum Likelihood Estimator and Leave-Out Estimator (within/across sections). 

The code is developed in relation to the paper Media-Polarisation in Denmark [here](https://github.com/Christianmoelby/Feature-extraction-Danish-news/blob/main/Media-Polarisation%20in%20Denmark.pdf). The paper contains appropriate citations of the components and on details of e.g. the dataset of faces to recognise. 

Using the processes below we have compiled a dataset of all relevant news articles and images from DR.dk and TV2.dk from 2015-2024, including all features listed above. Please reach out to get access to this dataset. 

## 🚀 Project Overview

This repository includes a series of Jupyter notebooks that together:

1. **Collect Danish news content**  
2. **Extract various features**—textual, visual, facial, and scene-based—from the content  
3. **Format feature representations** into token-based datasets for downstream analysis or modeling

X. **Estimate polarisation**

---

## 📂 Notebook Breakdown

Each notebook corresponds to one feature type or processing stage:

### 1. `Data generation.ipynb`  
- **Purpose**: Scrape and compile raw Danish news articles (text and images).  
- **Key steps**:  
  - Web scraping using Python (e.g., `requests`, `BeautifulSoup`)  
  - Storage of raw texts and images in structured folders/datasets  
- **Usage**: Run top-to-bottom to build raw data directory structure.

### 2. `Textual features.ipynb`  
- **Purpose**: Extract NLP-based features from article text.  
- **Techniques**:  
  - Tokenization  
  - mRAKE
- **Output**: Textual feature vectors per article ready for modeling.

### 3. `Facial_features.ipynb`  
- **Purpose**: Detect and vectorize faces present in images.  
- **Tools**: Deepface
- **Output**: Embedding vectors corresponding to detected faces in news images.

### 4. `Scene_tagging.ipynb`  
- **Purpose**: Analyze overall image scenes (e.g., indoor/outdoor, crowd, etc.).  
- **Tools**: OLLAMA
- **Output**: Scene classification tags as vectors for each image.

### 5. `Translate.ipynb`  
- **Purpose**: Fast translate of ENG image tags to DK.  
- **Output**: DK tokens.

### 6. `Compiling tokings.ipynb`  
- **Purpose**: Combine all previous feature outputs into a unified, tokenized dataset.  
- **Process**:  
  - Merge textual, visual, facial, and scene features  
  - Assign tokens/pointers to each feature group  
  - Output final structured dataset (CSV)

### 7. `MLE_Estimator.R, LO_Estimator_across_sections.R, LO_Estimator_within_sections.R`
- **Purpose**: Implement theoretical estimators of polarisation on token data generated from news.   
- **Output**:  
  - Unbaised polarisation measures across corpora or within sections, including randomly assigned control polarisation. 

---

## 📚 Recommended Citation

If you use this repository, please cite:

**Mølby, C. A., & Bremholm, K. A. H. (2025).** *Feature extraction – GitHub repository*. GitHub.  
https://github.com/Christianmoelby/Feature-extraction-Danish-news  

We reserve the right to continuously update and maintain the repository.

## ✉️ Contact: 
Christian Andersen Mølby: Christianmoelby@gmail.com

Knud Anton Højmark Bremholm: Knud.bremholm@gmail.com




