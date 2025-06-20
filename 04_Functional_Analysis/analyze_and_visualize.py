#analyze_and_visualize.py
import pandas as pd
import networkx as nx
import itertools
import os
import matplotlib.pyplot as plt
import re
import argparse

def analyze_and_visualize_modules(enrichment_folder: str, output_folder: str, cancer_type: str):
    """
    Processes all module enrichment files in a folder, groups them into functions,
    visualizes the networks, and saves a summary excel file.
    """
    
    print("--- step 3: functional grouping & visualization ---")

    if not os.path.isdir(enrichment_folder):
        print(f"error: enrichment folder '{enrichment_folder}' not found.")
        return
        
    if not os.path.exists(output_folder):
        print(f"Creating output folder: {output_folder}")
        os.makedirs(output_folder)

    enrichment_files = [f for f in os.listdir(enrichment_folder) if f.endswith(".csv")]
    if not enrichment_files:
        print("no enrichment .csv files found to analyze.")
        return
        
    #name summary file based on cancer type
    summary_excel_path = os.path.join(output_folder, f"{cancer_type}_Functions_Summary_Manual.xlsx")
    with pd.ExcelWriter(summary_excel_path, engine='xlsxwriter') as writer:
        #sort files numerically by module number
        sorted_files = sorted(enrichment_files, key=lambda x: int(re.search(r'_(\d+)_', x).group(1)))
        
        for file_name in sorted_files:
            match = re.search(r'_(\d+)_', file_name)
            if not match:
                continue
            module_num = match.group(1)
            module_name = f"Module {module_num}"
            print(f"--- analyzing {module_name} ---")

            file_path = os.path.join(enrichment_folder, file_name)
            pathways_df = pd.read_csv(file_path)
            
            significant_paths = pathways_df[
                (pathways_df['Fold_Enrichment'] > 1.5) & (pathways_df['FDR'] < 0.05)
            ].copy()

            if significant_paths.empty:
                print(f"No significant pathways to group for {module_name}.")
                pd.DataFrame(columns=["Functions", "Average FE"]).to_excel(
                    writer, sheet_name=module_name, index=False
                )
                continue

            print(f"Found {len(significant_paths)} significant pathways to analyze for {module_name}.")
            significant_paths['GeneSet'] = significant_paths['Genes'].apply(lambda x: set(str(x).split(', ')))

            pathway_to_geneset = pd.Series(significant_paths.GeneSet.values, index=significant_paths.Term).to_dict()

            pathway_graph = nx.Graph()
            for term, fe in pd.Series(significant_paths.Fold_Enrichment.values, index=significant_paths.Term).items():
                fe_value = fe if fe != float('inf') else significant_paths['Fold_Enrichment'][significant_paths['Fold_Enrichment'] != float('inf')].max() * 1.5
                pathway_graph.add_node(term, FE=fe_value)

            for path1, path2 in itertools.combinations(pathway_to_geneset.keys(), 2):
                gene_set1 = pathway_to_geneset[path1]
                gene_set2 = pathway_to_geneset[path2]
                
                num_common_genes = len(gene_set1.intersection(gene_set2))
                if num_common_genes > 0:
                    pathway_graph.add_edge(path1, path2, weight=num_common_genes)

            communities = list(nx.community.greedy_modularity_communities(pathway_graph, weight='weight'))
            sorted_communities = sorted([list(c) for c in communities], key=len, reverse=True)
            avg_fes = [sum(pathway_graph.nodes[p]['FE'] for p in c) / len(c) for c in sorted_communities]

            res_df = pd.DataFrame({'Functions': sorted_communities, 'Average FE': avg_fes})
            res_df.to_excel(writer, sheet_name=module_name, index=False)
            worksheet = writer.sheets[module_name]
            worksheet.set_column('A:A', 100)
            worksheet.set_column('B:B', 20)
            
            img_path = os.path.join(output_folder, f"{cancer_type}_Module_{module_num}_network.png")
            visualize_pathway_network(pathway_graph, sorted_communities, module_name, img_path)
            print(f"Saved visualization for {module_name} to {img_path}")

    print(f"--- step 3 complete. final summary saved to '{summary_excel_path}' ---")

def visualize_pathway_network(graph: nx.Graph, communities: list, title: str, output_image_file: str):
    if not communities or not graph.nodes():
        print(f"Skipping visualization for {title} as there are no communities or nodes.")
        return
        
    plt.figure(figsize=(24, 24))
    pos = nx.spring_layout(graph, k=0.9, iterations=60, seed=42)
    
    colors = plt.cm.get_cmap('tab20', len(communities))
    node_color_map = {}
    for i, comm in enumerate(communities):
        for node in comm:
            node_color_map[node] = colors(i)
            
    node_sizes = [graph.nodes[n]['FE'] * 25 for n in graph.nodes()]

    nx.draw_networkx_nodes(
        graph, pos,
        node_size=node_sizes,
        node_color=[node_color_map.get(node, 'lightgrey') for node in graph.nodes()]
    )
    nx.draw_networkx_edges(graph, pos, alpha=0.3, width=0.7)
    
    labels = {
        node: node.split('~')[-1].replace('_', ' ').title() 
        for node, size in zip(graph.nodes(), node_sizes) if size > 150
    }
    nx.draw_networkx_labels(graph, pos, labels=labels, font_size=10, font_weight='bold')
    
    plt.title(f"Functional Communities of Enriched Pathways for {title}", fontsize=28)
    plt.box(False)
    plt.tight_layout()
    plt.savefig(output_image_file, dpi=300, bbox_inches='tight')
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Visualize and group enriched pathways from modules.")
    parser.add_argument("--enrichment", required=True, help="Path to the folder with enrichment CSV files.")
    parser.add_argument("--output", required=True, help="Path to the output folder for saving images and summary.")
    parser.add_argument("--cancer", required=True, help="The name of the cancer type (e.g., 'BRCA') for file naming.")
    
    args = parser.parse_args()
        
    analyze_and_visualize_modules(
        enrichment_folder=args.enrichment,
        output_folder=args.output,
        cancer_type=args.cancer
    )