#run_all_cancers.py
import os
import subprocess
import pandas as pd
from collections import defaultdict

CANCER_TYPES = [
    "acc", "blca", "brca", "cesc", "coad", "esca", "hnsc", "kich",
    "kirc", "kirp", "laml", "lgg", "lihc", "luad", "lusc", "meso",
    "ov", "paad", "pcpg", "prad", "read", "sarc", "skcm", "stad",
    "tgct", "thca", 
    "thym", "ucec", "ucs", "uvm"
]

PATH_FIND_MODULES = r"02_Module_Discovery/find_modules.py"
PATH_RUN_ENRICHMENT = r"03_Pathway_Enrichment/run_enrichment.R"
PATH_ML_GROUPING = r"04_Functional_Analysis/ml_functional_grouping.py"
PATH_CREATE_SUMMARY = r"04_Functional_Analysis/create_final_summary.py"

NETWORK_INPUT_PATH = r"NetworkEdgelists"
BASE_INPUT_PATH = r"01_Input_Data"
BASE_MODULE_PATH = r"02_Module_Discovery"
BASE_ENRICHMENT_PATH = r"03_Pathway_Enrichment"
BASE_ANALYSIS_PATH = r"04_Functional_Analysis"
PAN_CANCER_OUTPUT_FOLDER = r"05_Pan_Cancer_Analysis" #new folder for final results

def run_command(command):
    """Helper function to run a command line process and check for errors."""
    print(f"\n>>> EXECUTING: {command}")
    try:
        subprocess.run(command, check=True, shell=True, capture_output=True, text=True)
        print(f">>> SUCCESS: Finished '{command.split()[0]}'")
    except subprocess.CalledProcessError as e:
        print(f"!!! ERROR: Command failed with exit code {e.returncode}")
        print(f"!!! FAILED COMMAND: {command}")
        print(f"--- STDOUT ---\n{e.stdout}")
        print(f"--- STDERR ---\n{e.stderr}")
        raise  #stop the pipeline if any step fails

def synthesize_pan_cancer_results():
    """
    Reads all the final summary tables and identifies conserved and
    cancer-specific miRNA-function relationships.
    """
    print("\n=========================================================")
    print("  STARTING PAN-CANCER SYNTHESIS")
    print("=========================================================\n")

    summary_files = [os.path.join(BASE_ANALYSIS_PATH, f) for f in os.listdir(BASE_ANALYSIS_PATH) if f.endswith("_Final_Paper_Table.csv")]
    
    #dictionaries to aggregate data across all cancers
    function_counts = defaultdict(int)
    function_to_mirnas = defaultdict(lambda: defaultdict(list))

    #read and aggregate data from each cancer's summary table
    for f_path in summary_files:
        cancer_type = os.path.basename(f_path).split('_')[0]
        df = pd.read_csv(f_path)
        for _, row in df.iterrows():
            #handle cases where functions might be "N/A" or "Not Significant"
            if pd.isna(row['Functions']) or "Significant" in row['Functions']:
                continue
            
            functions = [func.strip() for func in row['Functions'].split('+')]
            mirnas = [mirna.strip() for mirna in row['Top 3 miRNAs'].split(',\n')]
            
            for func in functions:
                function_counts[func] += 1
                function_to_mirnas[func][cancer_type] = mirnas

    #--- create pan cancer conserved functions tabl
    conserved_df = pd.DataFrame.from_dict(function_counts, orient='index', columns=['Cancer_Count'])
    conserved_df = conserved_df.sort_values(by='Cancer_Count', ascending=False).reset_index()
    conserved_df = conserved_df.rename(columns={'index': 'Functional_Group'})
    
    #filter for functions present in at least 5 cancers (example threshold)
    conserved_df_filtered = conserved_df[conserved_df['Cancer_Count'] >= 5]
    print("--- Top Conserved Pan-Cancer Functions ---")
    print(conserved_df_filtered.head(10))

    #--- create cancer-specific functions table ---
    specific_df = conserved_df[conserved_df['Cancer_Count'] == 1]
    #add the miRNA drivers for these specific functions
    specific_df['Cancer_Type'] = specific_df['Functional_Group'].apply(
        lambda func: list(function_to_mirnas[func].keys())[0]
    )
    specific_df['Top_3_miRNAs'] = specific_df['Functional_Group'].apply(
        lambda func: ", ".join(list(function_to_mirnas[func].values())[0])
    )
    print("\n--- Top Cancer-Specific Functions ---")
    print(specific_df.head(10))
    
    # --- save synthesis results to single .xlsx file ---
    # output_excel_path = os.path.join(PAN_CANCER_OUTPUT_FOLDER, "Pan_Cancer_Summary.xlsx")
    # with pd.ExcelWriter(output_excel_path, engine='xlsxwriter') as writer:
    #     conserved_df.to_excel(writer, sheet_name="Conserved_Functions", index=False)
    #     specific_df.to_excel(writer, sheet_name="Cancer_Specific_Functions", index=False)
    
    print(f"\n>>> Pan-cancer synthesis complete. Results saved to '{output_excel_path}'")

def main():
    """
    Main function to run the entire miRNA analysis pipeline
    for each specified cancer type.
    """
    #create the final output folder if it doesn't exist
    os.makedirs(PAN_CANCER_OUTPUT_FOLDER, exist_ok=True)
    
    #get the absolute path of the project's root directory
    project_root = os.path.abspath(os.path.dirname(__file__))

    for cancer_type in CANCER_TYPES:
        print(f"\n---------------------------------------------------------")
        print(f"  STARTING ANALYSIS FOR CANCER TYPE: {cancer_type.upper()}")
        print(f"---------------------------------------------------------\n")
        
        #helper function to create absolute paths
        def to_abs_path(rel_path):
            return os.path.abspath(os.path.join(project_root, rel_path)).replace('\\', '/')

        edge_list_file = to_abs_path(f"{NETWORK_INPUT_PATH}/{cancer_type}_EdgeList2.txt")
        gmt_folder = to_abs_path(f"{BASE_INPUT_PATH}/pathway_gmt_files")
        
        modules_file = to_abs_path(f"{BASE_MODULE_PATH}/{cancer_type}_Modules.txt")
        enrichment_folder = to_abs_path(f"{BASE_ENRICHMENT_PATH}/Enrichment_Results_{cancer_type}")
        
        ml_summary_folder = to_abs_path(f"{BASE_ANALYSIS_PATH}/Functional_Summary_ML_{cancer_type}")
        ml_summary_file = to_abs_path(f"{ml_summary_folder}/{cancer_type}_Functions_Summary_ML.xlsx")
        final_table_file = to_abs_path(f"{BASE_ANALYSIS_PATH}/{cancer_type}_Final_Paper_Table.csv")
        
        #create output directories
        os.makedirs(enrichment_folder.replace(project_root.replace('\\','/') + '/', ''), exist_ok=True)
        os.makedirs(ml_summary_folder.replace(project_root.replace('\\','/') + '/', ''), exist_ok=True)
        
        #--- execute the pipeline steps in order ---
        if not os.path.exists(edge_list_file):
            print(f"!!! WARNING: Edge list for {cancer_type} not found at {edge_list_file}. Skipping this cancer type.")
            continue

        #step 1: find_modules.py
        cmd_find_modules = f"python {PATH_FIND_MODULES} --edgelist \"{edge_list_file}\" --output \"{modules_file}\""
        run_command(cmd_find_modules)

        #step 2: run_enrichment.R
        cmd_run_enrichment = f"Rscript {PATH_RUN_ENRICHMENT} --modules \"{modules_file}\" --gmt \"{gmt_folder}\" --output \"{enrichment_folder}\""
        run_command(cmd_run_enrichment)

        #step 3: ml_functional_grouping.py
        cmd_ml_grouping = f"python {PATH_ML_GROUPING} --enrichment \"{enrichment_folder}\" --output \"{ml_summary_file}\""
        run_command(cmd_ml_grouping)
        
        #step 4: create_final_summary.py
        cmd_create_summary = f"python {PATH_CREATE_SUMMARY} --ml_summary \"{ml_summary_file}\" --modules \"{modules_file}\" --network \"{edge_list_file}\" --output \"{final_table_file}\""
        run_command(cmd_create_summary)

        print(f"\n--- COMPLETED ANALYSIS FOR {cancer_type.upper()} ---")

    #run final synthesis after all cancers are processed
    synthesize_pan_cancer_results()

if __name__ == "__main__":
    main()