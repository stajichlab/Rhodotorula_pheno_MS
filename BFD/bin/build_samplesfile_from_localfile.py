#!/usr/bin/env python3
"""
Build a samples CSV file for the BFD pipeline from local protein files.

Usage:
    module load taxonkit
    python bin/build_samplesfile_from_localfile.py -i input -o samples.csv

This script:
1. Scans input/pep/ for protein FASTA files
2. Extracts LOCUSTAG from the first protein header (>([^_]+)_ pattern)
3. Extracts species and strain from the filename
4. Uses taxonkit to look up NCBI taxon ID and taxonomy
5. Validates and fills missing taxonomy data from reference CSV
6. Generates a tab-delimited samples.csv file compatible with BFD pipeline

Output columns:
    ASMID: LOCAL_{filename_stem}
    SPECIES_IN: Full species name (genus + species from filename)
    STRAIN: Strain identifier (parsed from filename)
    NCBI_TAXONID: NCBI taxonomy ID (from reference CSV or taxonkit)
    BUSCO_LINEAGE: BUSCO lineage (from reference CSV)
    PHYLUM, CLASS, ORDER, FAMILY, etc: Full taxonomy (from reference CSV)
    LOCUSTAG: First component of protein header (pattern: >([^_]+)_)
    TRANSL_TABLE: Translation table (default: 1)

Example filename patterns:
    Rhodotorula_mucilaginosa_DBVPG_3045.proteins.fa
    -> Species: Rhodotorula mucilaginosa, Strain: DBVPG 3045

    Cystobasidium_sp._DVPG_10075.proteins.fa
    -> Species: Cystobasidium sp., Strain: DVPG 10075
"""

import os
import re
import sys
import subprocess
import argparse
from functools import lru_cache
from pathlib import Path
from typing import Optional, Dict, List, Tuple
import csv


def extract_locustag(fasta_file: Path) -> Optional[str]:
    """Extract LOCUSTAG from first protein header using pattern >([^_]+)_"""
    try:
        with open(fasta_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    match = re.match(r'>([^_]+)_', line)
                    if match:
                        return match.group(1)
                    else:
                        print(f"Warning: Could not extract LOCUSTAG from {fasta_file.name}: {line.strip()}",
                              file=sys.stderr)
                        return None
        print(f"Warning: No headers found in {fasta_file}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"Error reading {fasta_file}: {e}", file=sys.stderr)
        return None


def parse_filename(filename: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Extract species and strain from filename like:
    Rhodotorula_araucariae_NRRL_Y-17376.proteins.fa

    Returns: (species, strain) tuple
    """
    # Remove extension
    name = re.sub(r'\.(proteins|scaffolds)\.fa$', '', filename)

    # Try to parse genus_species_strain format
    parts = name.split('_')
    if len(parts) < 2:
        return None, None

    # "Genus_sp._clade_X_strain" -> species "Genus sp. clade X" (informal
    # clade designations used for undescribed Rhodotorula lineages), rest is
    # strain. Otherwise join first two parts as species (e.g. "Rhodotorula
    # araucariae").
    if parts[1] == 'sp.' and len(parts) > 3 and parts[2].lower() == 'clade':
        species = f"{parts[0]} sp. clade {parts[3]}"
        strain = ' '.join(parts[4:]) if len(parts) > 4 else None
    else:
        species = f"{parts[0]} {parts[1]}"
        strain = ' '.join(parts[2:]) if len(parts) > 2 else None

    return species, strain


@lru_cache(maxsize=None)
def lookup_taxonid(species: str, strain: Optional[str] = None) -> Tuple[Optional[int], Optional[str]]:
    """
    Use taxonkit to lookup NCBI taxon ID and full taxonomy for a species.
    Returns: (taxon_id, taxonomy_string) tuple

    Memoized on species (strain is unused in the lookup itself) since most
    samples share the same species and each call spawns two taxonkit
    subprocesses -- caching turns e.g. ~200 Rhodotorula mucilaginosa rows
    into a single taxonkit round-trip.
    """
    try:
        # taxonkit name2taxid reads names from stdin (a bare positional arg
        # is treated as a file path to read from, not a name), so the name
        # must be piped in via input=, not passed as an argv token.
        result = subprocess.run(
            ['taxonkit', 'name2taxid', '-s'],
            input=species + '\n',
            capture_output=True,
            text=True,
            timeout=10
        )

        if result.returncode == 0 and result.stdout.strip():
            # Parse output: species\ttaxid per line
            for line in result.stdout.strip().split('\n'):
                if line:
                    parts = line.split('\t')
                    if len(parts) >= 2 and parts[1]:
                        try:
                            taxid = int(parts[1])
                            # Get full taxonomy
                            lineage = get_lineage(taxid)
                            return taxid, lineage
                        except ValueError:
                            continue

        # Informal/undescribed names (e.g. "Rhodotorula sp. clade I") have no
        # NCBI entry of their own; fall back to the genus so NCBI_TAXONID
        # (required by the samples.csv schema) is still populated.
        genus = species.split()[0] if species else ''
        if genus and genus != species:
            return lookup_taxonid(genus)

        return None, None

    except subprocess.TimeoutExpired:
        return None, None
    except FileNotFoundError:
        return None, None
    except Exception:
        return None, None


def get_lineage(taxon_id: int) -> Optional[str]:
    """Get full taxonomic lineage for a taxon ID"""
    try:
        result = subprocess.run(
            ['taxonkit', 'lineage', str(taxon_id)],
            capture_output=True,
            text=True,
            timeout=10
        )

        if result.returncode == 0 and result.stdout.strip():
            # Output format: taxid\tlineage
            line = result.stdout.strip().split('\n')[0]
            parts = line.split('\t')
            if len(parts) >= 2:
                return parts[1]
    except Exception as e:
        print(f"Warning: Could not get lineage for taxon {taxon_id}: {e}", file=sys.stderr)

    return None


def parse_lineage_string(lineage: str) -> Dict[str, str]:
    """Parse taxonkit lineage string into taxonomy ranks"""
    ranks = {}

    if not lineage:
        return ranks

    # Lineage format: genus;species;kingdom;phylum;class;...
    # or similar - need to parse based on rank information
    # For now, return empty dict - this would need rank information from taxonkit

    return ranks


def load_existing_samples(csv_path: Path) -> Dict[str, Dict]:
    """Load existing samples CSV for reference/validation. Supports both tab and comma delimiters."""
    samples = {}
    try:
        with open(csv_path, 'r') as f:
            # Detect delimiter
            sample_line = f.readline()
            f.seek(0)
            delimiter = ',' if ',' in sample_line else '\t'

            reader = csv.DictReader(f, delimiter=delimiter)
            if reader.fieldnames:
                for row in reader:
                    if 'SPECIES_IN' in row and 'STRAIN' in row:
                        # Index by species (primary key)
                        species = row['SPECIES_IN'].strip()
                        strain = row['STRAIN'].strip()
                        key = f"{species}_{strain}" if strain else species
                        samples[key] = row
                        # Also index by just species for fuzzy matching
                        if species not in samples:
                            samples[species] = row
    except FileNotFoundError:
        pass
    except Exception as e:
        print(f"Warning: Could not load reference samples: {e}", file=sys.stderr)

    return samples


def process_protein_files(input_dir: Path, reference_csv: Optional[Path] = None,
                           busco_lineage: Optional[str] = None) -> List[Dict]:
    """
    Process all protein files in input directory and return sample records.
    Fills in missing taxonomy data from reference CSV when available.
    """
    pep_dir = input_dir / 'pep'

    if not pep_dir.exists():
        print(f"Error: {pep_dir} does not exist", file=sys.stderr)
        sys.exit(1)

    # Load reference samples for filling in taxonomy data
    reference_samples = {}
    if reference_csv and reference_csv.exists():
        reference_samples = load_existing_samples(reference_csv)

    # Case-insensitive species -> first-seen reference row, built once so
    # each of the 276 local files does an O(1) lookup instead of a fresh
    # linear scan over ~23k reference rows.
    species_index: Dict[str, Dict] = {}
    for row in reference_samples.values():
        key = row.get('SPECIES_IN', '').strip().lower()
        if key and key not in species_index:
            species_index[key] = row

    samples = []
    protein_files = sorted(pep_dir.glob('*.proteins.fa'))

    if not protein_files:
        print(f"Warning: No protein files found in {pep_dir}", file=sys.stderr)

    for prot_file in protein_files:
        print(f"Processing {prot_file.name}...", file=sys.stderr)

        # Extract LOCUSTAG
        locustag = extract_locustag(prot_file)
        if not locustag:
            print(f"Skipping {prot_file.name} - could not extract LOCUSTAG", file=sys.stderr)
            continue

        # Parse filename for species and strain
        species, strain = parse_filename(prot_file.name)
        if not species:
            print(f"Skipping {prot_file.name} - could not parse species/strain", file=sys.stderr)
            continue

        # Try to fill in taxonomy data from reference samples.
        # First try exact match with species and strain, then just species
        # (case-insensitive, since reference SPECIES_IN casing can vary).
        ref_key = f"{species}_{strain}" if strain else species
        ref = reference_samples.get(ref_key) or species_index.get(species.lower())

        # Only shell out to taxonkit (memoized per species) when the reference
        # CSV didn't already give us a taxon id -- the reference covers most
        # species in this dataset, so this keeps taxonkit calls rare.
        if ref and ref.get('NCBI_TAXONID'):
            taxid, lineage = None, None
        else:
            taxid, lineage = lookup_taxonid(species, strain)

        # Create sample record with basic info
        asmid = f"LOCAL_{prot_file.stem}"
        sample = {
            'ASMID': asmid,
            'SPECIES_IN': species,
            'STRAIN': strain or '',
            'BIOPROJECT': '',
            'NCBI_TAXONID': str(taxid) if taxid else '',
            'BUSCO_LINEAGE': '',
            'PHYLUM': '',
            'SUBPHYLUM': '',
            'CLASS': '',
            'SUBCLASS': '',
            'ORDER': '',
            'FAMILY': '',
            'GENUS': species.split()[0] if ' ' in species else '',
            'SPECIES': species,
            'TRANSL_TABLE': '1',
            'LOCUSTAG': locustag,
        }

        if ref:
            # Fill in missing values from reference
            for field in ['NCBI_TAXONID', 'BUSCO_LINEAGE', 'PHYLUM', 'SUBPHYLUM', 'CLASS', 'SUBCLASS', 'ORDER', 'FAMILY']:
                if not sample[field] and field in ref and ref[field]:
                    sample[field] = ref[field]

        # A caller-supplied lineage overrides whatever the reference CSV had
        # (e.g. forcing "basidiomycota" for a run that's all Basidiomycota).
        if busco_lineage:
            sample['BUSCO_LINEAGE'] = busco_lineage

        samples.append(sample)

    return samples


def write_samples_csv(samples: List[Dict], output_file: Path) -> None:
    """Write samples to CSV file (comma-delimited)"""
    if not samples:
        print("Warning: No samples to write", file=sys.stderr)
        return

    fieldnames = [
        'ASMID', 'SPECIES_IN', 'STRAIN', 'BIOPROJECT', 'NCBI_TAXONID',
        'BUSCO_LINEAGE', 'PHYLUM', 'SUBPHYLUM', 'CLASS', 'SUBCLASS',
        'ORDER', 'FAMILY', 'GENUS', 'SPECIES', 'TRANSL_TABLE', 'LOCUSTAG'
    ]

    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter=',', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(samples)

    print(f"Wrote {len(samples)} samples to {output_file}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description='Build BFD samples CSV from local protein files'
    )
    parser.add_argument(
        '-i', '--input',
        type=Path,
        default=Path.cwd(),
        help='Input directory containing pep/, dna/, cds/, gff3/ subdirectories (default: current dir)'
    )
    parser.add_argument(
        '-o', '--output',
        type=Path,
        default=None,
        help='Output samples CSV file (default: input_dir/samples.csv)'
    )
    parser.add_argument(
        '-r', '--reference',
        type=Path,
        default=Path('/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv'),
        help='Reference samples CSV for validation and filling missing values'
    )
    parser.add_argument(
        '-b', '--busco-lineage',
        type=str,
        default=None,
        help='Force BUSCO_LINEAGE to this value for every sample (overrides reference CSV)'
    )

    args = parser.parse_args()

    input_dir = args.input.resolve()
    if not input_dir.exists():
        print(f"Error: Input directory {input_dir} does not exist", file=sys.stderr)
        sys.exit(1)

    output_file = args.output or input_dir / 'samples.csv'

    # Process files
    samples = process_protein_files(input_dir, args.reference, args.busco_lineage)

    # Write output
    if samples:
        write_samples_csv(samples, output_file)
        print(f"Success: Created {output_file}")
    else:
        print("Error: No samples were generated", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
