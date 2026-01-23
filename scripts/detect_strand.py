#!/usr/bin/env python3

import sys
import subprocess
import argparse
import re

def parse_infer_experiment_output(output):
    """
    Parses the output of infer_experiment.py to determine strandedness.
    
    Expected output format from infer_experiment.py looks something like:
    This is PairEnd Data
    Fraction of reads explained by "1++,1--,2+-,2-+": 0.0156
    Fraction of reads explained by "1+-,1-+,2++,2--": 0.9844
    Unknown: 0.0000
    
    Or for SingleEnd:
    This is SingleEnd Data
    Fraction of reads explained by "++,--": 0.9844
    Fraction of reads explained by "+-,-+": 0.0156
    Unknown: 0.0000
    """
    
    lines = output.strip().split('\n')
    mode = None
    scores = {}
    
    for line in lines:
        if "PairEnd Data" in line:
            mode = "PE"
        elif "SingleEnd Data" in line:
            mode = "SE"
        
        # Parse fractions
        # PE: "1++,1--,2+-,2-+" -> Forward (FR) -> hisat2 --rna-strandness FR / F (if single)?? No wait.
        # Let's clarify HISAT2 vs RSeQC notation.
        
        # RSeQC:
        # 1++,1--,2+-,2-+  => read1 same strand as gene. This is "Forward" strandedness.
        # 1+-,1-+,2++,2--  => read1 opposite strand as gene. This is "Reverse" strandedness.
        
        # HISAT2 --rna-strandness:
        # F or R (Single End)
        # FR or RF (Paired End)
        
        # Mapping:
        # RSeQC "1++,1--,2+-,2-+" (Forward) => HISAT2 FR (or F for SE)
        # RSeQC "1+-,1-+,2++,2--" (Reverse) => HISAT2 RF (or R for SE)
        
        # SE:
        # "++,--" => Forward => HISAT2 F
        # "+-,-+" => Reverse => HISAT2 R
        
        if "Fraction of reads explained by" in line:
            parts = line.split(':')
            label = parts[0].strip().replace('Fraction of reads explained by ', '').replace('"', '')
            value = float(parts[1].strip())
            scores[label] = value

    if not mode:
        return None, "Unknown mode"

    # Define thresholds
    threshold = 0.6 # Lowered from 0.8
    
    detected_strand = None
    
    if mode == "PE":
        fwd_score = scores.get("1++,1--,2+-,2-+", 0.0)
        rev_score = scores.get("1+-,1-+,2++,2--", 0.0)
        
        # Primary check: meets absolute threshold
        if fwd_score >= threshold:
            detected_strand = "FR" # Forward
        elif rev_score >= threshold:
            detected_strand = "RF" # Reverse
        
        # Secondary check: Signal is strong (high ratio) even if below threshold (e.g. 0.6-0.7)
        # elif fwd_score > 0.6 and rev_score < 0.1:
        #      detected_strand = "FR"
        # elif rev_score > 0.6 and fwd_score < 0.1:
        #      detected_strand = "RF"

        elif fwd_score < 0.2 and rev_score < 0.2:
             # Likely unstranded if both are low? Or balanced? 
             # If balanced ~0.5 each, it's unstranded.
             # User asked to throw error if fails/ambiguous, or maybe they implied default unstranded?
             # User query: "throw an error and stop"
             pass
    
    elif mode == "SE":
        fwd_score = scores.get("++,--", 0.0)
        rev_score = scores.get("+-,-+", 0.0)
        
        if fwd_score >= threshold:
            detected_strand = "F" # Forward
        elif rev_score >= threshold:
            detected_strand = "R" # Reverse
        # Secondary check: Signal is strong (high ratio) even if below threshold (e.g. 0.6-0.7)
        # elif fwd_score > 0.6 and rev_score < 0.1:
        #      detected_strand = "F"
        # elif rev_score > 0.6 and fwd_score < 0.1:
        #      detected_strand = "R"
            
    if detected_strand:
        return detected_strand, None
    else:
        return None, f"Ambiguous result: {scores}"

def main():
    parser = argparse.ArgumentParser(description="Detect RNA strandedness using RSeQC")
    parser.add_argument("--bam", required=True, help="Input BAM file")
    parser.add_argument("--bed", required=True, help="Reference BED file")
    
    args = parser.parse_args()
    
    # Run infer_experiment.py
    cmd = ["infer_experiment.py", "-r", args.bed, "-i", args.bam]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        detected_strand, error = parse_infer_experiment_output(result.stdout)
        
        if detected_strand:
            print(detected_strand)
            sys.exit(0)
        else:
            print(f"Error: Could not detect strandedness. {error}", file=sys.stderr)
            # print stdout for debugging
            print(result.stdout, file=sys.stderr)
            sys.exit(1)
            
    except subprocess.CalledProcessError as e:
        print(f"Error running infer_experiment.py: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
