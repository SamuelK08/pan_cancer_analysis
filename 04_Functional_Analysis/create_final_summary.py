#create_final_summary.py
import os
import pandas as pd
import networkx as nx
import re
import argparse 

def create_final_summary_table(
    ml_summary_path: str,
    modules_path: str,
    network_path: str,
    output_csv_path: str
):
    """
    Combines the ML-generated functions with the top miRNAs for each module
    to create a final, paper-ready summary table.
    """
    print("--- Step 4: Creating the Final Summary Table ---")

    # --- 1. load all necessary data files ---
    print(f"Loading ML summary from: {ml_summary_path}")
    try:
        if not os.path.exists(ml_summary_path) or os.path.getsize(ml_summary_path) < 100:
            print("ML summary file is empty or missing. Creating an empty summary table.")
            #create an empty CSV to signify completion and prevent a crash
            pd.DataFrame(columns=["Community", "Functions", "Top 3 miRNAs"]).to_csv(output_csv_path, index=False)
            return 0 #exit
        
        ml_summary_sheets = pd.read_excel(ml_summary_path, sheet_name=None)
    except Exception as e:
        print(f"ERROR: Could not read Excel file: {e}")
        print("Creating an empty summary table to allow pipeline to continue.")
        pd.DataFrame(columns=["Community", "Functions", "Top 3 miRNAs"]).to_csv(output_csv_path, index=False)
        return 1


    print(f"Loading modules from: {modules_path}")
    try:
        with open(modules_path, 'r') as f:
            module_lines = f.readlines()
        if not module_lines:
            print("Warning: Modules file is empty. Cannot generate summary.")
            pd.DataFrame(columns=["Community", "Functions", "Top 3 miRNAs"]).to_csv(output_csv_path, index=False)
            return 0
        module_list = [line.strip().split(', ') for line in module_lines]
    except FileNotFoundError:
        print(f"ERROR: Modules file not found at '{modules_path}'")
        return 1

    print(f"Loading network from: {network_path}")
    try:
        G = nx.read_edgelist(network_path)
    except FileNotFoundError:
        print(f"ERROR: Network edge list not found at '{network_path}'")
        return 1

    final_summary_data = []

    #--- 2. iterate through each module to assemble the final table ---
    num_modules = len(module_list)
    print(f"Processing {num_modules} modules...")

    for i in range(num_modules):
        module_num = i + 1
        module_name = f"Module {module_num}"

        #--- 2a. get the summarized functions for this module ---
        functions_string = "Not Significant" 
        if module_name in ml_summary_sheets:
            module_df = ml_summary_sheets[module_name]
            if not module_df.empty:
                top_functions = module_df['Functional_Group'].head(3).tolist()
                functions_string = " + ".join(top_functions)
        
        if not functions_string or functions_string.strip() == "":
            functions_string = "N/A"

        #--- 2b. calculate the top 3 miRNAs based on local degree ---
        top_mirnas_string = "N/A"
        module_nodes = module_list[i]
        
        module_subgraph = G.subgraph(module_nodes)
        mirna_nodes = [node for node in module_nodes if node.startswith('hsa-')]

        if mirna_nodes:
            mirna_degrees = dict(module_subgraph.degree(mirna_nodes))
            if mirna_degrees: #check if dictionary is not empty
                sorted_mirnas = sorted(mirna_degrees.items(), key=lambda item: item[1], reverse=True)
                top_3_mirnas = [mirna[0] for mirna in sorted_mirnas[:3]]
                top_mirnas_string = ",\n".join(top_3_mirnas)

        #--- 2c. append the row to our summary list ---
        final_summary_data.append({
            "Community": module_num,
            "Functions": functions_string,
            "Top 3 miRNAs": top_mirnas_string
        })

    #--- 3. create and save the final dataframe ---
    final_df = pd.DataFrame(final_summary_data)
    final_df.to_csv(output_csv_path, index=False)
    
    print("\n--- Final Summary Table ---")
    print(final_df.to_string()) #use to_string() to ensure all rows are printed
    print(f"\nSuccessfully saved the final summary to '{output_csv_path}'")
    print("--- Project summary complete! ---")
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create a final summary table for a cancer type.")
    
    parser.add_argument("--ml_summary", required=True, help="Path to the ML-generated functions Excel file.")
    parser.add_argument("--modules", required=True, help="Path to the modules definition text file.")
    parser.add_argument("--network", required=True, help="Path to the original network edge list file.")
    parser.add_argument("--output", required=True, help="Path to save the final output CSV summary table.")

    args = parser.parse_args()

    #run the function and exit with a code
    exit_code = create_final_summary_table(
        ml_summary_path=args.ml_summary,
        modules_path=args.modules,
        network_path=args.network,
        output_csv_path=args.output
    )
    if exit_code != 0:
        exit(exit_code)