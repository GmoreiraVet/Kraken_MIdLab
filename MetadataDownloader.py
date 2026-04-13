from Bio import Entrez
from Bio import SeqIO
import json
import time
import os

# =========================
# CONFIGURATION
# =========================
Entrez.email = "your_email@example.com"

output_dir = "/home/viroicbas2023/Documents/Gmoreira/EnteroVirusCarraçasUropigeaCustom/Enterovirus/PhyloTree/WholeGenome/30RandomFullgenomes/GB_Files"

download_format = "gb"   # Options: "gb", "fasta", "xml", etc.
save_sequences = True    # Save downloaded records to files
delay = 0.4              # NCBI rate limiting safety

# =========================
# INPUT ACCESSIONS
# =========================
accessions = [
    "AJQ21577", "QBP14846", "AYP73841", "XSC61406", "BAD89965",
    "AHN91979", "QSV20191", "ACH87535", "UXK63315", "ALB72985",
    "ABC69262", "WGJ63450", "AJD77366", "AJD77360", "WIL60417",
    "WIL60428", "WVR18621", "AHX39363", "BBD18130", "BCU46238",
    "WHA31286", "QJA10586", "AEK98146", "AEK98147", "ABC69251",
    "ABC69255", "QED91176", "QED91177", "AUI41040", "BAN79276",
    "AVQ54966", "AGO64489", "AGO64482", "XYL45182", "XYL45181",
    "ABC69250", "ABC69249", "AVQ54963", "AFD03543", "AFD03545",
    "XBS34000", "XBS33999", "AVQ54965", "AVQ55056", "QIN85535",
    "ABE60737", "AVQ54964"
]

# =========================
# SETUP
# =========================
os.makedirs(output_dir, exist_ok=True)

# Mapping Entrez rettype → SeqIO format + extension
FORMAT_MAP = {
    "gb": ("genbank", ".gb"),
    "gbwithparts": ("genbank", ".gb"),
    "fasta": ("fasta", ".fasta")
}

seqio_format, file_ext = FORMAT_MAP.get(download_format, (None, f".{download_format}"))

all_metadata = []

# =========================
# MAIN LOOP
# =========================
for acc in accessions:
    try:
        print(f"Fetching {acc}...")

        handle = Entrez.efetch(
            db="protein",
            id=acc,
            rettype=download_format,
            retmode="text"
        )

        record = None

        # -------------------------
        # CASE 1: SeqIO-supported formats
        # -------------------------
        if seqio_format:
            try:
                record = SeqIO.read(handle, seqio_format)
            finally:
                handle.close()

            if save_sequences and record:
                out_file = os.path.join(output_dir, f"{acc}{file_ext}")
                SeqIO.write(record, out_file, seqio_format)

        # -------------------------
        # CASE 2: Raw formats (XML, ASN.1, etc.)
        # -------------------------
        else:
            data = handle.read()
            handle.close()

            if save_sequences:
                out_file = os.path.join(output_dir, f"{acc}{file_ext}")
                with open(out_file, "w") as f:
                    f.write(data)

        # -------------------------
        # METADATA EXTRACTION (GenBank only)
        # -------------------------
        if seqio_format == "genbank" and record:

            # --- Source qualifiers ---
            source_quals = {}
            for feature in record.features:
                if feature.type == "source":
                    for key in [
                        "host", "country", "collection_date",
                        "isolate", "strain", "serotype", "note"
                    ]:
                        if key in feature.qualifiers:
                            source_quals[key] = feature.qualifiers[key][0]

            # --- Publications ---
            publications = []
            for ref in record.annotations.get("references", []):
                pub = {
                    "title": ref.title or None,
                    "authors": ref.authors or None,
                    "journal": ref.journal or None,
                    "pubmed": ref.pubmed_id.strip() if ref.pubmed_id else None,
                }

                if pub["title"] and pub["title"] != "Direct Submission":
                    publications.append(pub)

            # --- Metadata entry ---
            entry = {
                "accession": acc,
                "description": record.description,
                "organism": record.annotations.get("organism"),
                "taxonomy": record.annotations.get("taxonomy", []),
                "host": source_quals.get("host"),
                "country": source_quals.get("country"),
                "collection_date": source_quals.get("collection_date"),
                "isolate": source_quals.get("isolate"),
                "strain": source_quals.get("strain"),
                "serotype": source_quals.get("serotype"),
                "note": source_quals.get("note"),
                "publications": publications,
            }

            all_metadata.append(entry)

        print(f"  OK")
        time.sleep(delay)

    except Exception as e:
        print(f"  ERROR fetching {acc}: {e}")

# =========================
# SAVE METADATA
# =========================
if all_metadata:
    metadata_path = os.path.join(output_dir, "metadata.json")
    with open(metadata_path, "w") as f:
        json.dump(all_metadata, f, indent=2, default=str)

    print(f"\nMetadata saved to {metadata_path}")

print("\nDone.")