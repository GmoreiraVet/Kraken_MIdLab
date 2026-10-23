
#!/usr/bin/env python3
"""
Fetch collection dates and source metadata from GenBank for every
accession in a FASTA file.

Requirements:
    pip install biopython

Run in VS Code:
    1. Edit the CONFIG section below.
    2. Run this script.

Outputs:
    genbank_metadata.tsv
    dates.tsv

The dates.tsv file contains:
    accession.version<TAB>collection_date

Note:
    FASTA tip names should match accession.version for downstream
    phylogenetic analyses. Accessions that cannot be identified
    are reported as MISSING_ACCESSION.
"""

# =============================================================================
# CONFIG: EDIT THESE SETTINGS
# =============================================================================

CONFIG = {
    # Input FASTA file
    "fasta_file": (
        "/home/viroicbas2023/Documents/Gmoreira/Carraças/Carraças_Metagenomica_Gui/VOLZSHKOESEPARADO/Filogenia/M_Segment/M_Segment_Database_FInal_ConsensusPT.fasta"
    ),

    # Output directory. Use "." for the current working directory.
    "output_dir": "/home/viroicbas2023/Documents/Gmoreira/Carraças/Carraças_Metagenomica_Gui/VOLZSHKOESEPARADO/MolecularClock_Analysis/M_Segment",

    # NCBI Entrez email (required by NCBI)
    "email": "gmoreiravet@gmail.com",

    # Optional NCBI API key. Use "" if you do not have one.
    "api_key": "",

    # Request handling
    "max_retries": 3,
    "retry_delay_seconds": 2,
    "request_delay_seconds": 0.4,

    # Output filenames
    "metadata_filename": "genbank_metadata.tsv",
    "dates_filename": "dates.tsv",

    # FASTA identifier parsing
    # Supports common nucleotide accession formats, including:
    # NC_000001.1, OP123456.1, ABC12345.1
    "accession_pattern": (
        r"^((?:[A-Z]{1,2}_\d{5,9}|"
        r"[A-Z]{1,4}\d{5,9}|"
        r"[A-Z]{2}\d{6,8}|"
        r"[A-Z]{4}\d{8,10})"
        r"(?:\.\d+)?)"
    ),
}

# =============================================================================
# IMPORTS
# =============================================================================

import re
import sys
import time
from pathlib import Path

from Bio import Entrez, SeqIO


# =============================================================================
# INITIALISATION
# =============================================================================

def validate_config():
    """Validate configuration before starting network requests."""

    fasta_path = Path(CONFIG["fasta_file"]).expanduser()

    if not fasta_path.is_file():
        raise FileNotFoundError(
            f"Input FASTA file not found:\n{fasta_path}\n"
            "Edit CONFIG['fasta_file'] at the top of the script."
        )

    if not CONFIG["email"] or CONFIG["email"] == "you@example.org":
        raise ValueError(
            "Set your email address in CONFIG['email'] before running."
        )

    output_dir = Path(CONFIG["output_dir"]).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)

    return fasta_path, output_dir


def get_accession(record_id, accession_regex):
    """Extract a GenBank accession.version from a FASTA record ID."""

    match = accession_regex.match(record_id)
    return match.group(1) if match else None


# =============================================================================
# GENBANK RETRIEVAL
# =============================================================================

def fetch_genbank_record(accession):
    """Fetch a GenBank record, retrying failed requests."""

    for attempt in range(1, CONFIG["max_retries"] + 1):
        handle = None

        try:
            handle = Entrez.efetch(
                db="nuccore",
                id=accession,
                rettype="gb",
                retmode="text",
            )

            record = SeqIO.read(handle, "genbank")
            return record

        except Exception as error:
            print(
                f"  Attempt {attempt}/{CONFIG['max_retries']} failed "
                f"for {accession}: {error}",
                file=sys.stderr,
            )

            if attempt < CONFIG["max_retries"]:
                time.sleep(CONFIG["retry_delay_seconds"])

        finally:
            if handle is not None:
                handle.close()

    return None


def extract_source_metadata(record):
    """Extract metadata from the GenBank source feature."""

    source_feature = next(
        (
            feature
            for feature in record.features
            if feature.type == "source"
        ),
        None,
    )

    if source_feature is None:
        qualifiers = {}
    else:
        qualifiers = source_feature.qualifiers

    def get_qualifier(key):
        values = qualifiers.get(key, [])
        return values[0] if values else ""

    collection_date = get_qualifier("collection_date")
    country = get_qualifier("country")
    host = get_qualifier("host")
    strain_isolate = (
        get_qualifier("strain") or get_qualifier("isolate")
    )

    status = "OK" if collection_date else "NO_DATE"

    return (
        collection_date,
        country,
        host,
        strain_isolate,
        status,
    )


# =============================================================================
# MAIN WORKFLOW
# =============================================================================

def main():
    fasta_path, output_dir = validate_config()

    Entrez.email = CONFIG["email"]

    if CONFIG["api_key"]:
        Entrez.api_key = CONFIG["api_key"]

    accession_regex = re.compile(CONFIG["accession_pattern"])

    metadata_path = output_dir / CONFIG["metadata_filename"]
    dates_path = output_dir / CONFIG["dates_filename"]

    print(f"Input FASTA: {fasta_path}")
    print(f"Output directory: {output_dir}")
    print("Reading FASTA records...")

    records = list(SeqIO.parse(str(fasta_path), "fasta"))

    if not records:
        raise ValueError("The FASTA file contains no sequences.")

    rows = []
    total = len(records)

    print(f"Found {total} sequences.\n")

    for index, record in enumerate(records, start=1):
        name = record.id
        accession = get_accession(name, accession_regex)

        print(f"[{index}/{total}] {name}")

        if accession is None:
            print("  No accession detected.")
            rows.append(
                (
                    name,
                    "",
                    "",
                    "",
                    "",
                    "",
                    "MISSING_ACCESSION",
                )
            )
            continue

        gb_record = fetch_genbank_record(accession)

        if gb_record is None:
            rows.append(
                (
                    name,
                    accession,
                    "",
                    "",
                    "",
                    "",
                    "FETCH_FAILED",
                )
            )
            continue

        (
            collection_date,
            country,
            host,
            strain_isolate,
            status,
        ) = extract_source_metadata(gb_record)

        rows.append(
            (
                name,
                accession,
                collection_date,
                country,
                host,
                strain_isolate,
                status,
            )
        )

        print(
            f"  Accession: {accession} | "
            f"Date: {collection_date or 'NOT AVAILABLE'} | "
            f"Status: {status}"
        )

        # Respect NCBI request-rate limits.
        time.sleep(CONFIG["request_delay_seconds"])

    # -------------------------------------------------------------------------
    # Write complete metadata table
    # -------------------------------------------------------------------------

    columns = (
        "name",
        "accession",
        "collection_date",
        "country",
        "host",
        "strain_isolate",
        "status",
    )

    with metadata_path.open("w", encoding="utf-8", newline="") as output:
        output.write("\t".join(columns) + "\n")

        for row in rows:
            output.write("\t".join(row) + "\n")

    # -------------------------------------------------------------------------
    # Write date file for phylogenetic tools
    # -------------------------------------------------------------------------

    with dates_path.open("w", encoding="utf-8", newline="") as output:
        for row in rows:
            name, accession, date = row[0], row[1], row[2]

            if accession and date:
                output.write(f"{accession}\t{date}\n")

    # -------------------------------------------------------------------------
    # Summary
    # -------------------------------------------------------------------------

    successful = sum(row[6] == "OK" for row in rows)
    no_date = sum(row[6] == "NO_DATE" for row in rows)
    failed = sum(row[6] == "FETCH_FAILED" for row in rows)
    missing = sum(row[6] == "MISSING_ACCESSION" for row in rows)

    print("\n" + "=" * 60)
    print("COMPLETED")
    print("=" * 60)
    print(f"Total sequences:              {total}")
    print(f"Records with collection date: {successful}")
    print(f"Records without date:         {no_date}")
    print(f"Failed GenBank requests:      {failed}")
    print(f"Missing accessions:           {missing}")
    print(f"\nMetadata table: {metadata_path}")
    print(f"Date file:      {dates_path}")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, FileNotFoundError) as error:
        sys.exit(f"\nERROR: {error}")
