#!/usr/bin/env python3

import pandas as pd
import sys
import os
import matplotlib.pyplot as plt

def load_data(file_path):
    try:
        # Try reading as TSV
        return pd.read_csv(file_path, sep='\t')
    except:
        # Try CSV if TSV fails
        return pd.read_csv(file_path)

def clean_v_call(v_call):
    if pd.isna(v_call):
        return "None"
    # Take the first call if comma separated, remove allele (*01)
    return v_call.split(',')[0].split('*')[0]

def compare_results(pipeline_file, reference_file, outdir):
    print(f"Loading pipeline result: {pipeline_file}")
    df_pipe = load_data(pipeline_file)
    
    print(f"Loading reference data: {reference_file}")
    df_ref = load_data(reference_file)
    
    # Ensure sequence_id is the key
    # If cell_id is the key, rename or ensure column exists
    key = 'sequence_id'
    if key not in df_pipe.columns:
        print(f"Error: {key} not in pipeline output")
        # Check for other likely keys
        if 'sequence_id' in df_pipe.columns: key = 'sequence_id'
        elif 'query_id' in df_pipe.columns: key = 'query_id' # Immcantation sometimes uses query_id
    
    if key not in df_pipe.columns:
         print(f"FATAL: Could not find key column in pipeline output. Columns: {df_pipe.columns}")
         sys.exit(1)
         
    # Rename pipeline key to match ref if needed (assuming Ref uses sequence_id)
    # Ref: sequence_id
    
    # Merge
    merged = pd.merge(df_pipe, df_ref, on=key, how='inner', suffixes=('_pipe', '_ref'))
    
    total_shared = len(merged)
    print(f"Shared sequences: {total_shared}")
    
    # Compare V-calls (Gene level)
    # Pipeline might produce "IGHV1-18*01", Ref "IGHV1-18*01"
    # We strip allele for robustness
    merged['v_call_pipe_clean'] = merged['v_call_pipe'].apply(clean_v_call)
    merged['v_call_ref_clean'] = merged['v_call_ref'].apply(clean_v_call)
    
    merged['v_match'] = merged['v_call_pipe_clean'] == merged['v_call_ref_clean']
    v_accuracy = merged['v_match'].mean()
    
    # Compare J-calls
    merged['j_call_pipe_clean'] = merged['j_call_pipe'].apply(clean_v_call)
    merged['j_call_ref_clean'] = merged['j_call_ref'].apply(clean_v_call)
    merged['j_match'] = merged['j_call_pipe_clean'] == merged['j_call_ref_clean']
    j_accuracy = merged['j_match'].mean()

    # Compare CDR3 AA / Junction AA
    # Pipeline output might have junction_aa instead of cdr3_aa
    # Reference has cdr3_aa AND junction_aa
    
    col_pipe = 'cdr3_aa'
    if 'cdr3_aa' not in merged.columns and 'cdr3_aa_pipe' not in merged.columns:
        if 'junction_aa' in df_pipe.columns:
             col_pipe = 'junction_aa'
        elif 'junction_aa_pipe' in merged.columns:
             col_pipe = 'junction_aa' # handled by suffix logic if name was same?
    
    # Actually, after merge we have _pipe and _ref suffixes if collision
    # If both had 'cdr3_aa', we get cdr3_aa_pipe and cdr3_aa_ref
    # If Pipe has 'junction_aa' and Ref has 'junction_aa', we get junction_aa_pipe ...
    
    # Let's try to find the columns dynamically
    pipe_cdr3 = merged['cdr3_aa_pipe'] if 'cdr3_aa_pipe' in merged.columns else (merged['junction_aa_pipe'] if 'junction_aa_pipe' in merged.columns else None)
    ref_cdr3 = merged['cdr3_aa_ref'] if 'cdr3_aa_ref' in merged.columns else (merged['cdr3_aa'] if 'cdr3_aa' in merged.columns else None)
    
    # If Ref file was loaded and NOT renamed during merge (because Pipe didn't have same col), it keeps original name.
    # Ref has 'cdr3_aa'. Pipe has 'junction_aa'. They don't collide.
    # So merged has 'junction_aa' and 'cdr3_aa'.
    
    if pipe_cdr3 is None:
        # Fallback: maybe it didn't collide?
        if 'junction_aa' in merged.columns: pipe_cdr3 = merged['junction_aa']
        
    if ref_cdr3 is None:
        if 'cdr3_aa' in merged.columns: ref_cdr3 = merged['cdr3_aa']
        
    if pipe_cdr3 is None or ref_cdr3 is None:
        print("Warning: Could not find CDR3/Junction AA columns for comparison.")
        cdr3_accuracy = 0.0
    else:
        pipe_cdr3 = pipe_cdr3.fillna("")
        ref_cdr3 = ref_cdr3.fillna("")
        merged['cdr3_match'] = pipe_cdr3 == ref_cdr3
        cdr3_accuracy = merged['cdr3_match'].mean()
    
    print(f"V-Gene Accuracy: {v_accuracy:.2%}")
    print(f"J-Gene Accuracy: {j_accuracy:.2%}")
    print(f"CDR3 Accuracy: {cdr3_accuracy:.2%}")

    # Save summary
    summary = pd.DataFrame({
        'Metric': ['Shared_Sequences', 'V_Gene_Accuracy', 'J_Gene_Accuracy', 'CDR3_Accuracy'],
        'Value': [total_shared, v_accuracy, j_accuracy, cdr3_accuracy]
    })
    summary.to_csv(f"{outdir}/comparison_summary.csv", index=False)
    
    # Plot
    metrics = ['V_Gene', 'J_Gene', 'CDR3']
    values = [v_accuracy, j_accuracy, cdr3_accuracy]
    
    plt.figure(figsize=(8, 6))
    plt.bar(metrics, values, color=['skyblue', 'lightgreen', 'salmon'])
    plt.ylim(0, 1.05)
    plt.title(f"Comparison with Reference (n={total_shared})")
    plt.ylabel("Accuracy (Exact Match)")
    
    for i, v in enumerate(values):
        plt.text(i, v + 0.02, f"{v:.1%}", ha='center')
        
    plt.savefig(f"{outdir}/comparison_plot.png")

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: compare_results.py <pipeline_csv> <reference_csv> <outdir>")
        sys.exit(1)
        
    compare_results(sys.argv[1], sys.argv[2], sys.argv[3])
