#find_modules.py
import networkx as nx
from cdlib import algorithms
import os
import argparse 

def find_and_save_modules(edge_list_path: str, output_path: str):
    print("--- step 1: module discovery (using cdlib/louvain algorithm) ---")
    if not os.path.exists(edge_list_path):
        print(f"ERROR: Input edge list '{edge_list_path}' not found.")
        return 1 #return an error

    output_dir = os.path.dirname(output_path)
    if not os.path.exists(output_dir):
        print(f"Creating output directory: {output_dir}")
        os.makedirs(output_dir, exist_ok=True)

    print(f"Loading graph from: {edge_list_path}")
    G = nx.read_edgelist(edge_list_path)
    if not G.nodes():
        print("ERROR: The graph is empty. Please check your edge list file.")
        return 1 #return an error
        
    largest_cc = max(nx.connected_components(G), key=len)
    G_lcc = G.subgraph(largest_cc)
    print(f"Graph loaded. LCC has {G_lcc.number_of_nodes()} nodes and {G_lcc.number_of_edges()} edges.")

    print("Running louvain algorithm to find modules...")
    try:
        coms = algorithms.louvain(G_lcc, randomize=False) #add randomize=False for deterministic results
        mod_list = coms.communities
        
        if not mod_list:
            print("WARNING: Louvain algorithm did not find any communities.")
            #still create an empty file so the pipeline doesn't break
            with open(output_path, 'w') as f:
                pass #create empty file
            return 0

    except Exception as e:
        print(f"ERROR: Louvain algorithm failed: {e}")
        return 1

    print(f"Louvain identified {len(mod_list)} modules.")

    #write the modules to the output file
    with open(output_path, 'w') as f:
        for module_nodes in mod_list:
            nodes = sorted(module_nodes)
            line = ", ".join(map(str, nodes))
            f.write(line + "\n")
    
    print(f"Successfully saved modules to '{output_path}'")
    print("--- step 1 complete ---")
    return 0 #return success code

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find network modules from an edge list.")
    parser.add_argument("--edgelist", required=True, help="Path to the input edge list file.")
    parser.add_argument("--output", required=True, help="Path to save the output modules file.")
    
    args = parser.parse_args()

    #exit with a non-zero code if the function fails
    exit_code = find_and_save_modules(edge_list_path=args.edgelist, output_path=args.output)
    if exit_code != 0:
        exit(exit_code)