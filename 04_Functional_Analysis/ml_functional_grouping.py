#ml_functional_grouping.py
import os
import pandas as pd
import numpy as np
import re
from sentence_transformers import SentenceTransformer
from sklearn.cluster import KMeans
from sklearn.feature_extraction.text import TfidfVectorizer
import argparse 

def clean_term_name(term):
    if '~' in term:
        term = term.split('~', 1)[1]
    term = re.sub(r'^(GO|KEGG|REACTOME|WP|BIOCARTA|PID)_', '', term, flags=re.IGNORECASE)
    term = term.replace('_', ' ').lower()
    return term

def get_cluster_name(terms_in_cluster, vectorizer, feature_names):
    cluster_text = " ".join(terms_in_cluster)
    if not cluster_text.strip():
        return "Unnamed Cluster"
    tfidf_matrix = vectorizer.transform([cluster_text])
    non_zero_indices = tfidf_matrix.toarray()[0].nonzero()[0]
    if len(non_zero_indices) == 0:
        return "Broadly Associated Terms"
    sorted_indices = sorted(non_zero_indices, key=lambda i: tfidf_matrix[0, i], reverse=True)
    top_indices = sorted_indices[:3]
    top_words = [feature_names[i] for i in top_indices]
    return ", ".join(word.capitalize() for word in top_words)

def group_functions_with_ml(enrichment_folder: str, output_path: str):
    print("--- step 3 (ML): automated functional grouping ---")

    if not os.path.isdir(enrichment_folder):
        print(f"Error: Enrichment folder '{enrichment_folder}' not found.")
        return 1
        
    enrichment_files = [f for f in os.listdir(enrichment_folder) if f.endswith(".csv")]
    if not enrichment_files:
        print("No enrichment .csv files found.")
        #create an empty excel file to signify completion
        pd.DataFrame().to_excel(output_path)
        return 0

    print("Loading NLP model (this may take a moment)...")
    model = SentenceTransformer('all-MiniLM-L6-v2')

    with pd.ExcelWriter(output_path, engine='xlsxwriter') as writer:
        sorted_files = sorted(enrichment_files, key=lambda x: int(re.search(r'_(\d+)_', x).group(1)))
        
        for file_name in sorted_files:
            match = re.search(r'_(\d+)_', file_name)
            if not match: continue
            
            module_num = match.group(1)
            module_name = f"Module {module_num}"
            print(f"--- Analyzing {module_name} with ML ---")

            file_path = os.path.join(enrichment_folder, file_name)
            pathways_df = pd.read_csv(file_path)

            significant_paths = pathways_df[
                (pathways_df['Fold_Enrichment'] > 2.0) & (pathways_df['FDR'] < 0.05)
            ].copy()
            
            if np.isinf(significant_paths['Fold_Enrichment']).any():
                max_fe = significant_paths.loc[np.isfinite(significant_paths['Fold_Enrichment']), 'Fold_Enrichment'].max()
                if pd.isna(max_fe): max_fe = 100 # handle case where all are Inf
                significant_paths['Fold_Enrichment'].replace(np.inf, max_fe * 1.5, inplace=True)

            if len(significant_paths) < 10: 
                print(f"Skipping {module_name}, not enough significant pathways ({len(significant_paths)} found).")
                pd.DataFrame().to_excel(writer, sheet_name=module_name)
                continue

            significant_paths['Clean_Term'] = significant_paths['Term'].apply(clean_term_name)
            pathway_names = significant_paths['Clean_Term'].tolist()

            print(f"Generating embeddings for {len(pathway_names)} pathways...")
            embeddings = model.encode(pathway_names, show_progress_bar=False)

            num_pathways = len(significant_paths)
            k = int(np.sqrt(num_pathways / 2)) #adjusted heuristic for better grouping
            k = max(2, min(k, 15)) 
            if k > len(significant_paths): k = len(significant_paths) #k cannot be > num samples
            
            print(f"Clustering into {k} functional groups...")
            kmeans = KMeans(n_clusters=k, random_state=42, n_init='auto')
            significant_paths['Cluster'] = kmeans.fit_predict(embeddings)

            print("Generating automatic names for clusters...")
            vectorizer = TfidfVectorizer(stop_words='english', max_df=0.9, min_df=1)
            vectorizer.fit(pathway_names)
            feature_names = vectorizer.get_feature_names_out()
            
            cluster_names = {}
            for i in range(k):
                terms_in_cluster = significant_paths[significant_paths['Cluster'] == i]['Clean_Term'].tolist()
                cluster_names[i] = get_cluster_name(terms_in_cluster, vectorizer, feature_names)

            significant_paths['Functional_Group'] = significant_paths['Cluster'].map(cluster_names)

            summary = significant_paths.groupby('Functional_Group').agg(
                Pathway_Count=('Term', 'size'),
                Avg_FE=('Fold_Enrichment', 'mean'),
                Example_Pathways=('Clean_Term', lambda x: " | ".join(x.head(3).str.title()))
            ).reset_index().sort_values(by='Avg_FE', ascending=False)
            
            summary.to_excel(writer, sheet_name=module_name, index=False)
            worksheet = writer.sheets[module_name]
            worksheet.set_column('A:A', 40)
            worksheet.set_column('B:B', 15)
            worksheet.set_column('C:C', 15)
            worksheet.set_column('D:D', 100)

    print(f"--- ML analysis complete. Final summary saved to '{output_path}' ---")
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Group enriched pathways into functional themes using ML.")
    parser.add_argument("--enrichment", required=True, help="Path to the folder with enrichment CSV files.")
    parser.add_argument("--output", required=True, help="Full path for the output Excel summary file.")
    
    args = parser.parse_args()
        
    exit_code = group_functions_with_ml(
        enrichment_folder=args.enrichment,
        output_path=args.output
    )
    if exit_code != 0:
        exit(exit_code)