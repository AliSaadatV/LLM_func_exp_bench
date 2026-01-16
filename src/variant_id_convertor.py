"""
Variant annotation via VariantValidator (rest.variantvalidator.org) with rsID fallback via Ensembl VEP.

Returns for BOTH hg38 and hg37:
  - HGVSc (transcript-level)
  - HGVSp (predicted protein consequence strings; may be multiple)
  - rsIDs (from VV when present; otherwise resolved via VEP colocated_variants)
  - chr_pos_ref_alt (built from VV VCF fields)

Inputs supported:
  - chrom/pos/ref/alt
  - or a single `variant` string (HGVS c/p/g, "17-50198002-C-A", "GRCh38:17:...", "chr17:50198002C>A",
    or "chr17_50198002_C_A", etc.)
"""

from __future__ import annotations

import re
from typing import Any, Dict, List, Optional, Tuple
from urllib.parse import quote

import requests

# ----------------------------
# VariantValidator config
# ----------------------------

VV_BASE = "https://rest.variantvalidator.org"
VV_PATH = "/VariantValidator/variantvalidator/{genome_build}/{variant_description}/{select_transcripts}"

# Transcript selection modes supported by VariantValidator (commonly used in their UI/docs)
ALLOWED_SELECT_TRANSCRIPTS = {
    "mane_select",
    "mane",
    "select",
    "refseq_select",
    "all",
    "raw",
}

# ----------------------------
# Ensembl VEP config (rsID fallback)
# ----------------------------

ENSEMBL = {
    "hg38": "https://rest.ensembl.org",
    "hg37": "https://grch37.rest.ensembl.org",
}
VEP_REGION_ENDPOINT = "/vep/homo_sapiens/region"

# ----------------------------
# Parsers / regex
# ----------------------------

# Accept common chr_pos_ref_alt-ish patterns like:
#   chr17_50198002_C_A, 17:50198002:C:A, 17-50198002-C-A, etc.
_CHRPOS_RE = re.compile(
    r"^(?:chr)?(?P<chrom>[0-9]{1,2}|X|Y|M|MT)[\s:_-]+(?P<pos>\d+)[\s:_-]+(?P<ref>[ACGTN-]+)[\s:_-]+(?P<alt>[ACGTN-]+)$",
    re.IGNORECASE,
)
_RS_RE = re.compile(r"\brs(\d+)\b", re.IGNORECASE)


# ----------------------------
# Small utilities
# ----------------------------

def _dedup_keep_order(xs: List[str]) -> List[str]:
    seen = set()
    out: List[str] = []
    for x in xs:
        if x and x not in seen:
            seen.add(x)
            out.append(x)
    return out


def _normalize_rsid(s: str) -> Optional[str]:
    m = _RS_RE.search(s or "")
    if not m:
        return None
    return f"rs{m.group(1)}"


def _normalize_variant_description(
    variant: Optional[str],
    chrom: Optional[str],
    pos: Optional[int],
    ref: Optional[str],
    alt: Optional[str],
) -> str:
    """
    Normalize input into something VariantValidator accepts.

    Priority:
      1) explicit chrom/pos/ref/alt -> "CHR-POS-REF-ALT"
      2) HGVS-ish strings -> pass through
      3) VV pseudo-VCF-ish strings -> pass through
      4) chr_pos_ref_alt -> convert to "CHR-POS-REF-ALT"
      5) otherwise pass through unchanged
    """
    if chrom is not None and pos is not None and ref is not None and alt is not None:
        c = str(chrom).strip()
        c = c[3:] if c.lower().startswith("chr") else c
        c = c.upper().replace("MT", "M")
        return f"{c}-{int(pos)}-{str(ref).upper()}-{str(alt).upper()}"

    if not variant:
        raise ValueError("Either `variant` or all of (chrom, pos, ref, alt) must be provided.")

    v = variant.strip()

    # Likely HGVS (transcript/genomic/protein) -> pass through
    if any(token in v for token in (":c.", ":p.", ":g.", ":m.", ":n.", " c.", " p.", " g.")):
        return v

    # Common VV pseudo-VCF formats -> pass through
    # - 17-50198002-C-A
    # - 17:50198002:C:A
    # - GRCh38:17:50198002:C:A
    # - chr17:50198002C>A
    if (v.count("-") == 3) or (v.count(":") in (3, 4)) or (">" in v):
        return v

    # chr_pos_ref_alt -> convert
    m = _CHRPOS_RE.match(v)
    if m:
        c = m.group("chrom").upper().replace("MT", "M")
        return f"{c}-{m.group('pos')}-{m.group('ref').upper()}-{m.group('alt').upper()}"

    return v


def _vv_get(
    genome_build: str,
    variant_description: str,
    select_transcripts: str,
    timeout: int,
) -> Dict[str, Any]:
    url = VV_BASE + VV_PATH.format(
        genome_build=quote(genome_build, safe=""),
        variant_description=quote(variant_description, safe=""),
        select_transcripts=quote(select_transcripts, safe=""),
    )
    r = requests.get(url, params={"content-type": "application/json"}, timeout=timeout)
    r.raise_for_status()
    return r.json()


def _find_main_record(payload: Dict[str, Any]) -> Tuple[str, Dict[str, Any]]:
    """
    VV response usually looks like:
      {
        "<variant_key>": {...main...},
        "flag": "...",
        "metadata": {...}
      }
    """
    for k, v in payload.items():
        if k in ("flag", "metadata"):
            continue
        if isinstance(v, dict):
            return k, v
    raise ValueError("Could not locate main VariantValidator record in response.")


def _extract_hgvsp(main: Dict[str, Any]) -> List[str]:
    """
    VV often provides hgvs_predicted_protein_consequence as a dict or string.
    """
    out: List[str] = []
    pc = main.get("hgvs_predicted_protein_consequence")
    if isinstance(pc, dict):
        for v in pc.values():
            if isinstance(v, str) and v:
                out.append(v)
    elif isinstance(pc, str) and pc:
        out.append(pc)
    return _dedup_keep_order(out)


def _extract_rsids_anywhere(obj: Any) -> List[str]:
    """
    VV does not reliably include rsIDs. When it does, scan recursively for rs####.
    """
    found: List[str] = []

    def walk(x: Any) -> None:
        if isinstance(x, dict):
            for kk, vv in x.items():
                if isinstance(kk, str):
                    rs = _normalize_rsid(kk)
                    if rs:
                        found.append(rs)
                walk(vv)
        elif isinstance(x, list):
            for it in x:
                walk(it)
        elif isinstance(x, str):
            rs = _normalize_rsid(x)
            if rs:
                found.append(rs)

    walk(obj)
    return _dedup_keep_order(found)


def _vcf_to_chr_pos_ref_alt(vcf: Dict[str, Any]) -> List[str]:
    """
    Convert VV VCF dict {chr,pos,ref,alt} to ["chrX_pos_ref_alt", ...].
    ALT may be comma-separated.
    """
    chrom = str(vcf.get("chr", "")).strip()
    pos = str(vcf.get("pos", "")).strip()
    ref = str(vcf.get("ref", "")).strip()
    alt = str(vcf.get("alt", "")).strip()

    if not (chrom and pos and ref and alt):
        return []

    if not chrom.lower().startswith("chr"):
        chrom = "chr" + chrom
    # normalize to "chr" + whatever follows
    chrom = "chr" + chrom[3:]

    alts = [a.strip() for a in alt.split(",") if a.strip()]
    return [f"{chrom}_{pos}_{ref}_{a}" for a in alts]


def _extract_build_locus(main: Dict[str, Any], wanted: str) -> Dict[str, Any]:
    """
    wanted: 'hg38' or 'hg37'
    VV primary_assembly_loci keys can include: hg38/grch38 and hg19/grch37.
    """
    pal = main.get("primary_assembly_loci", {})
    if not isinstance(pal, dict):
        return {"key": None, "hgvs_genomic": None, "vcf": None, "chr_pos_ref_alt": []}

    if wanted == "hg38":
        keys = ["hg38", "grch38"]
    elif wanted == "hg37":
        keys = ["hg19", "grch37"]  # hg19 corresponds to GRCh37
    else:
        raise ValueError("wanted must be 'hg38' or 'hg37'")

    for k in keys:
        block = pal.get(k)
        if isinstance(block, dict):
            hgvs_g = block.get("hgvs_genomic_description")
            vcf = block.get("vcf") if isinstance(block.get("vcf"), dict) else None
            chrpos = _vcf_to_chr_pos_ref_alt(vcf or {})
            return {
                "key": k,
                "hgvs_genomic": hgvs_g if isinstance(hgvs_g, str) else None,
                "vcf": vcf,
                "chr_pos_ref_alt": chrpos,
            }

    return {"key": None, "hgvs_genomic": None, "vcf": None, "chr_pos_ref_alt": []}


# ----------------------------
# VEP rsID fallback
# ----------------------------

def _vep_rsids_for_variant(build: str, chrom: str, pos: int, ref: str, alt: str, timeout: int) -> List[str]:
    """
    Ask Ensembl VEP (region endpoint) for colocated rsIDs.

    Uses a whitespace-delimited VCF-like line:
      CHROM POS ID REF ALT QUAL FILTER INFO
    """
    server = ENSEMBL[build]
    variant_str = f"{chrom} {pos} . {ref} {alt} . . ."

    payload = {
        "variants": [variant_str],
        # We only need colocated_variants; keep it simple.
        "hgvs": 0,
        "pick": 0,
    }
    headers = {"Content-Type": "application/json", "Accept": "application/json"}

    r = requests.post(server + VEP_REGION_ENDPOINT, headers=headers, json=payload, timeout=timeout)
    r.raise_for_status()
    data = r.json()
    if not data:
        return []

    entry = data[0]
    rsids: List[str] = []

    for cv in entry.get("colocated_variants", []) or []:
        vid = cv.get("id")
        if isinstance(vid, str):
            rs = _normalize_rsid(vid)
            if rs:
                rsids.append(rs)
        for x in (cv.get("ids") or []):
            if isinstance(x, str):
                rs = _normalize_rsid(x)
                if rs:
                    rsids.append(rs)

    return _dedup_keep_order(rsids)


def _resolve_rsids_from_vv_vcf_blocks(vv_hg38_block: Dict[str, Any], vv_hg37_block: Dict[str, Any], timeout: int) -> Dict[str, List[str]]:
    """
    Use VV-provided VCF coordinates to query VEP and recover rsIDs.
    Returns {"hg38": [...], "hg37": [...]}.
    """
    out = {"hg38": [], "hg37": []}

    for build, block in (("hg38", vv_hg38_block), ("hg37", vv_hg37_block)):
        vcf = (block or {}).get("vcf") or {}
        if not isinstance(vcf, dict):
            continue

        chrom = vcf.get("chr")
        pos = vcf.get("pos")
        ref = vcf.get("ref")
        alt = vcf.get("alt")

        if not (chrom and pos and ref and alt):
            continue

        chrom_s = str(chrom).strip()
        if chrom_s.lower().startswith("chr"):
            chrom_s = chrom_s[3:]

        # ALT can be comma-separated; query each ALT
        alts = [a.strip() for a in str(alt).split(",") if a.strip()]
        rs_for_build: List[str] = []
        for a in alts:
            rs_for_build.extend(
                _vep_rsids_for_variant(
                    build=build,
                    chrom=chrom_s,
                    pos=int(pos),
                    ref=str(ref),
                    alt=a,
                    timeout=timeout,
                )
            )
        out[build] = _dedup_keep_order(rs_for_build)

    return out


# ----------------------------
# Public API
# ----------------------------

def vv_annotate_variant(
    chrom: Optional[str] = None,
    pos: Optional[int] = None,
    ref: Optional[str] = None,
    alt: Optional[str] = None,
    variant: Optional[str] = None,
    genome_build: str = "GRCh38",
    select_transcripts: str = "mane_select",
    resolve_rsids: bool = True,
    timeout: int = 30,
) -> Optional[Dict[str, Any]]:
    """
    VariantValidator-based annotation with optional rsID fallback via Ensembl VEP.

    Returns a dict with:
      - rsids: union of rsids (VV + fallback)
      - rsids_source: {"vv": [...], "vep": {"hg38":[...], "hg37":[...]}}
      - hgvsc
      - hgvsp_all
      - gene_symbol
      - transcript_accession
      - hg38: {key, hgvs_genomic, vcf, chr_pos_ref_alt}
      - hg37: {key, hgvs_genomic, vcf, chr_pos_ref_alt}
      - warnings
      - flag, metadata
    """
    if select_transcripts not in ALLOWED_SELECT_TRANSCRIPTS:
        raise ValueError(f"select_transcripts must be one of: {sorted(ALLOWED_SELECT_TRANSCRIPTS)}")

    vd = _normalize_variant_description(variant=variant, chrom=chrom, pos=pos, ref=ref, alt=alt)

    payload = _vv_get(
        genome_build=genome_build,
        variant_description=vd,
        select_transcripts=select_transcripts,
        timeout=timeout,
    )
    if not payload:
        return None

    _, main = _find_main_record(payload)

    hgvsc = main.get("hgvs_transcript_variant") if isinstance(main.get("hgvs_transcript_variant"), str) else None
    hgvsp_all = _extract_hgvsp(main)
    gene_symbol = main.get("gene_symbol") if isinstance(main.get("gene_symbol"), str) else None
    warnings = main.get("validation_warnings") if isinstance(main.get("validation_warnings"), list) else []

    transcript_accession = None
    if hgvsc and ":" in hgvsc:
        transcript_accession = hgvsc.split(":", 1)[0]

    # rsIDs from VV (often empty)
    rsids_vv = _extract_rsids_anywhere(payload)

    hg38_block = _extract_build_locus(main, "hg38")
    hg37_block = _extract_build_locus(main, "hg37")

    rsids_vep_by_build = {"hg38": [], "hg37": []}
    rsids_final = rsids_vv[:]

    # Fallback: resolve via VEP colocated_variants using VV’s per-build VCF coords
    if resolve_rsids and not rsids_final:
        try:
            rsids_vep_by_build = _resolve_rsids_from_vv_vcf_blocks(
                vv_hg38_block=hg38_block,
                vv_hg37_block=hg37_block,
                timeout=timeout,
            )
            rsids_final = _dedup_keep_order(rsids_vep_by_build["hg38"] + rsids_vep_by_build["hg37"])
        except requests.RequestException:
            # Don't fail the whole annotation if VEP is down; just return without rsIDs
            rsids_vep_by_build = {"hg38": [], "hg37": []}

    return {
        "input": {
            "variant_description": vd,
            "genome_build_called": genome_build,
            "select_transcripts": select_transcripts,
        },
        "rsids": rsids_final,
        "rsids_source": {
            "vv": rsids_vv,
            "vep": rsids_vep_by_build,
        },
        "hgvsc": hgvsc,
        "hgvsp_all": hgvsp_all,
        "gene_symbol": gene_symbol,
        "transcript_accession": transcript_accession,
        "hg38": hg38_block,
        "hg37": hg37_block,
        "warnings": warnings,
        "flag": payload.get("flag"),
        "metadata": payload.get("metadata"),
    }


# ----------------------------
# Example usage
# ----------------------------
# if __name__ == "__main__":
#     examples = [
#         # chrom/pos/ref/alt
#         dict(chrom="17", pos=50198002, ref="C", alt="A"),
#         # HGVS
#         dict(variant="NM_000088.3:c.589G>T"),
#         # pseudo-VCF / flexible forms
#         dict(variant="17-50198002-C-A"),
#         dict(variant="GRCh38:17:50198002:C:A"),
#         dict(variant="chr17:50198002C>A"),
#         dict(variant="chr17_50198002_C_A"),
#     ]

#     for ex in examples:
#         res = vv_annotate_variant(**ex, select_transcripts="mane_select", resolve_rsids=True, timeout=30)
#         print("\nINPUT:", ex)
#         print("rsids:", res["rsids"] if res else None)
#         if res:
#             print("hg38 chr_pos_ref_alt:", res["hg38"]["chr_pos_ref_alt"])
#             print("hg37 chr_pos_ref_alt:", res["hg37"]["chr_pos_ref_alt"])
#             print("hgvsc:", res["hgvsc"])
#             print("hgvsp_all:", res["hgvsp_all"])
