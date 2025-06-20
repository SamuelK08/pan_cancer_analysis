#synthesize_pan_cancer.py
import os
import pandas as pd
from collections import defaultdict

def synthesize_pan_cancer_results(analysis_folder: str, output_excel_path: str):
    """
    Reads all final summary tables from each cancer analysis and synthesizes
    the results to find conserved and cancer-specific miRNA-function relationships.
    """
    print("\n=========================================================")
    print("  STARTING PAN-CANCER SYNTHESIS")
    print("=========================================================\n")

    # --- 1. find and load all individual cancer summary tables ---
    try:
        summary_files = [os.path.join(analysis_folder, f) for f in os.listdir(analysis_folder) if f.endswith("_Final_Paper_Table.csv")]
        if not summary_files:
            print(f"ERROR: No '*_Final_Paper_Table.csv' files found in '{analysis_folder}'.")
            return
        print(f"Found {len(summary_files)} cancer summary tables to analyze.")
    except FileNotFoundError:
        print(f"ERROR: Analysis folder not found at '{analysis_folder}'")
        return

    # --- 2. aggregate data across all cancer types ---
    function_to_cancer_mirnas = defaultdict(dict)

    for f_path in summary_files:
        cancer_type = os.path.basename(f_path).split('_')[0].upper()
        try:
            df = pd.read_csv(f_path)
            #handle cases where the CSV might be empty
            if df.empty:
                continue
        except pd.errors.EmptyDataError:
            print(f"Warning: Skipping empty summary file for {cancer_type}.")
            continue

        for _, row in df.iterrows():
            #skip rows with no significant function
            if pd.isna(row['Functions']) or "Significant" in row['Functions'] or "N/A" in row['Functions']:
                continue
            
            #split function groups like "A + B" into individual functions "A" and "B"
            functions = [func.strip() for func in row['Functions'].split('+')]
            
            # clean up miRNA list
            if pd.isna(row['Top 3 miRNAs']): continue
            mirnas = [mirna.strip() for mirna in row['Top 3 miRNAs'].split(',\n')]
            
            for func in functions:
                function_to_cancer_mirnas[func][cancer_type] = mirnas

    #--- 3. analyze aggregated data to find conserved functions ---
    print("Analyzing for conserved functions...")
    conserved_data = []
    for func, cancer_mirnas in function_to_cancer_mirnas.items():
        cancer_count = len(cancer_mirnas)
        #get a flat list of all miRNAs associated with this function across all cancers
        all_mirnas = [mirna for mirna_list in cancer_mirnas.values() for mirna in mirna_list]
        #find the top 3 most common miRNAs driving this function across all cancers
        top_pan_cancer_mirnas = pd.Series(all_mirnas).value_counts().nlargest(3).index.tolist()
        
        conserved_data.append({
            "Functional_Group": func,
            "Cancer_Count": cancer_count,
            "Top_Pan_Cancer_miRNAs": ", ".join(top_pan_cancer_mirnas),
            "Present_In_Cancers": ", ".join(sorted(cancer_mirnas.keys()))
        })

    conserved_df = pd.DataFrame(conserved_data).sort_values(by="Cancer_Count", ascending=False)

    #--- 4. analyze aggregated data to find cancer-specific functions ---
    print("Analyzing for cancer-specific functions...")
    specific_data = []
    for func, cancer_mirnas in function_to_cancer_mirnas.items():
        if len(cancer_mirnas) == 1: #found in only one cancer type
            cancer_type = list(cancer_mirnas.keys())[0]
            top_mirnas = list(cancer_mirnas.values())[0]
            
            specific_data.append({
                "Cancer_Type": cancer_type,
                "Functional_Group": func,
                "Top_3_miRNAs": ", ".join(top_mirnas)
            })
            
    specific_df = pd.DataFrame(specific_data).sort_values(by=["Cancer_Type", "Functional_Group"])

    #--- 5. save the final synthesis to a single Excel file ---
    print(f"Saving pan-cancer summary to: {output_excel_path}")
    with pd.ExcelWriter(output_excel_path, engine='xlsxwriter') as writer:
        conserved_df.to_excel(writer, sheet_name="Pan-Cancer_Conserved_Functions", index=False)
        specific_df.to_excel(writer, sheet_name="Cancer_Specific_Functions", index=False)
        
        # Auto-adjust column widths for readability
        for sheet_name in writer.sheets:
            worksheet = writer.sheets[sheet_name]
            for idx, col in enumerate(conserved_df if "Conserved" in sheet_name else specific_df):
                series = (conserved_df if "Conserved" in sheet_name else specific_df)[col]
                max_len = max((series.astype(str).map(len).max(), len(str(series.name)))) + 2
                worksheet.set_column(idx, idx, max_len)

    print("\n=========================================================")
    print("  PAN-CANCER SYNTHESIS COMPLETE!")
    print("=========================================================\n")


if __name__ == "__main__":
    ANALYSIS_FOLDER = r"04_Functional_Analysis"
    
    PAN_CANCER_FOLDER = r"05_Pan_Cancer_Analysis"
    os.makedirs(PAN_CANCER_FOLDER, exist_ok=True)
    OUTPUT_EXCEL_FILE = os.path.join(PAN_CANCER_FOLDER, "Pan_Cancer_miRNA_Function_Summary.xlsx")

    #run the synthesis function
    synthesize_pan_cancer_results(
        analysis_folder=ANALYSIS_FOLDER,
        output_excel_path=OUTPUT_EXCEL_FILE
    )