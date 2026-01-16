from typing import List, Dict, Any, Optional, Literal, Sequence
import re

def format_variant_for_prompt(variant_data: Dict[str, Any]) -> str:
    """
    Format the output of vv_annotate_variant into a string for use in prompts.
    
    Args:
        variant_data: Dictionary output from vv_annotate_variant
        
    Returns:
        Formatted string containing:
        - gene_symbol
        - vcf style (hg38: chr_pos_ref_alt)
        - vcf style (hg37: chr_pos_ref_alt)
        - genomics (hg38:g. ...)
        - genomics (hg37:g. ...)
        - rsids
        - hgvsc (just the part after :)
        - hgvsp (just the part after :)
    """
    parts = []
    
    # Gene symbol
    gene_symbol = variant_data.get('gene_symbol', 'N/A')
    parts.append(f"Gene: {gene_symbol}")
    
    # hg38 chr_pos_ref_alt (no "VCF" label)
    hg38_vcf = variant_data.get('hg38', {}).get('chr_pos_ref_alt', [])
    hg38_vcf_str = hg38_vcf[0] if hg38_vcf else 'N/A'
    parts.append(f"chr_pos_ref_alt (hg38): {hg38_vcf_str}")
    
    # hg37 chr_pos_ref_alt (no "VCF" label)
    hg37_vcf = variant_data.get('hg37', {}).get('chr_pos_ref_alt', [])
    hg37_vcf_str = hg37_vcf[0] if hg37_vcf else 'N/A'
    parts.append(f"chr_pos_ref_alt (hg37): {hg37_vcf_str}")
    
    # Genomics hg38 (extract part after 'g.')
    hg38_genomic = variant_data.get('hg38', {}).get('hgvs_genomic', '')
    if 'g.' in hg38_genomic:
        hg38_genomic_part = hg38_genomic.split('g.')[1]
    else:
        hg38_genomic_part = 'N/A'
    parts.append(f"HGVSg (hg38): g.{hg38_genomic_part}")
    
    # Genomics hg37 (extract part after 'g.')
    hg37_genomic = variant_data.get('hg37', {}).get('hgvs_genomic', '')
    if 'g.' in hg37_genomic:
        hg37_genomic_part = hg37_genomic.split('g.')[1]
    else:
        hg37_genomic_part = 'N/A'
    parts.append(f"HGVSg (hg37): g.{hg37_genomic_part}")
    
    # RSIDs
    rsids = variant_data.get('rsids', [])
    rsids_str = ', '.join(rsids) if rsids else 'N/A'
    parts.append(f"rsID: {rsids_str}")
    
    # HGVSc (extract part after ':')
    hgvsc = variant_data.get('hgvsc', '')
    hgvsc_part = hgvsc.split(':', 1)[1] if ':' in hgvsc else hgvsc if hgvsc else 'N/A'
    parts.append(f"HGVSc: {hgvsc_part}")
    
    # HGVSp (extract part after ':' from each entry, keep one 3-letter and one 1-letter instance)
    hgvsp_all = variant_data.get('hgvsp_all', [])
    if hgvsp_all:
        hgvsp_parts = []
        for hgvsp_entry in hgvsp_all:
            if ':' in hgvsp_entry:
                hgvsp_parts.append(hgvsp_entry.split(':', 1)[1])
            else:
                hgvsp_parts.append(hgvsp_entry)
        
        # Filter to entries with 3-letter amino acids (pattern like p.(Gly197Cys) or similar)
        # 3-letter amino acid codes: Ala, Arg, Asn, Asp, Cys, Gln, Glu, Gly, His, Ile, Leu, Lys, Met, Phe, Pro, Ser, Thr, Trp, Tyr, Val
        three_letter_pattern = re.compile(r'\([A-Z][a-z]{2}\d+[A-Z][a-z]{2}\)')
        three_letter_hgvsp = [hgvsp for hgvsp in hgvsp_parts if three_letter_pattern.search(hgvsp)]
        
        # Filter to entries with 1-letter amino acids (pattern like p.(G197C) or similar)
        # 1-letter amino acid codes: A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V
        one_letter_pattern = re.compile(r'\([A-Z]\d+[A-Z]\)')
        one_letter_hgvsp = [hgvsp for hgvsp in hgvsp_parts if one_letter_pattern.search(hgvsp)]
        
        # Get one instance of each (first one found)
        three_letter_str = three_letter_hgvsp[0] if three_letter_hgvsp else 'N/A'
        one_letter_str = one_letter_hgvsp[0] if one_letter_hgvsp else 'N/A'
        
        # Remove parentheses from both
        three_letter_str = three_letter_str.replace('(', '').replace(')', '') if three_letter_str != 'N/A' else 'N/A'
        one_letter_str = one_letter_str.replace('(', '').replace(')', '') if one_letter_str != 'N/A' else 'N/A'
        
        # Add both to output
        parts.append(f"HGVSp (3-letter): {three_letter_str}")
        parts.append(f"HGVSp (1-letter): {one_letter_str}")
    else:
        parts.append(f"HGVSp (3-letter): N/A")
        parts.append(f"HGVSp (1-letter): N/A")
    
    return ', '.join(parts)