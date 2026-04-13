from Bio import Entrez
import json
import time
import os

# --- CONFIGURE THESE ---
Entrez.email = "your_email@example.com"
output_dir = "/home/viroicbas2023/Downloads"  # <-- Add your path here
# -----------------------

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

os.makedirs(output_dir, exist_ok=True)

all_metadata = {}

for acc in accessions:
    try:
        print(f"Fetching metadata for {acc}...")

        # Step 1: get UID
        search_handle = Entrez.esearch(db="protein", term=acc)
        search_record = Entrez.read(search_handle)
        search_handle.close()

        uid = search_record["IdList"][0] if search_record["IdList"] else None
        if not uid:
            print(f"  No UID found for {acc}, skipping.")
            continue

        # Step 2: esummary metadata
        summary_handle = Entrez.esummary(db="protein", id=uid)
        summary = Entrez.read(summary_handle)
        summary_handle.close()
        summary_data = dict(summary[0])

        # Step 3: efetch GenBank record to extract host, country, collection_date
        fetch_handle = Entrez.efetch(db="protein", id=uid, rettype="gb", retmode="text")
        gb_text = fetch_handle.read()
        fetch_handle.close()

        # Parse key fields from the flat file text
        host, country, collection_date, isolate = None, None, None, None
        for line in gb_text.splitlines():
            line = line.strip()
            if line.startswith("/host="):
                host = line.split("=", 1)[1].strip('"')
            elif line.startswith("/country="):
                country = line.split("=", 1)[1].strip('"')
            elif line.startswith("/collection_date="):
                collection_date = line.split("=", 1)[1].strip('"')
            elif line.startswith("/isolate="):
                isolate = line.split("=", 1)[1].strip('"')

        all_metadata[acc] = {
            "accession":       acc,
            "title":           summary_data.get("Title", ""),
            "organism":        summary_data.get("Organism", ""),
            "taxid":           str(summary_data.get("TaxId", "")),
            "length":          str(summary_data.get("Length", "")),
            "create_date":     str(summary_data.get("CreateDate", "")),
            "update_date":     str(summary_data.get("UpdateDate", "")),
            "source_db":       summary_data.get("SourceDb", ""),
            "accession_ver":   summary_data.get("AccessionVersion", ""),
            "host":            host,
            "country":         country,
            "collection_date": collection_date,
            "isolate":         isolate,
        }

        print(f"  OK: {summary_data.get('Title', 'N/A')}")
        time.sleep(0.4)

    except Exception as e:
        print(f"  ERROR fetching {acc}: {e}")

# Save everything to one JSON file
out_path = os.path.join(output_dir, "metadata.json")
with open(out_path, "w") as f:
    json.dump(all_metadata, f, indent=2, default=str)

print(f"\nAll metadata saved to {out_path}")