"""
Messy Lab, step 2: reference lists.

Single source for the compound, cell-line and alias tables. Running it writes the CSVs next to it,
which pipeline/build.py loads into the database. Edit the lists here, then run:

    python3 pipeline/reference/make_reference.py

Checked 27 Sep 2026:
  * PubChem CIDs from the PubChem name lookup.
  * Cellosaurus IDs and diseases from the Cellosaurus API.
  * Every alias marked as seen in a study's GEO files was confirmed to appear in that study's raw
    files by check_reference.py. Aliases seen only in a paper say so in seen_in_field.
"""
import csv
from pathlib import Path

HERE = Path(__file__).parent

# compound_id, preferred_name, compound_class, is_vehicle, pubchem_cid, notes
COMPOUNDS = [
    ("OHT", "4-hydroxytamoxifen", "SERM", False, "449459",
     "Active metabolite of tamoxifen. Used by S1, S3, S4, S5. S1 specifies the trans isomer."),
    ("TAM_UNSPECIFIED", "tamoxifen (form not stated)", "SERM", False, "",
     "S2 says only 'tamoxifen (Sigma-Aldrich)'. May be the parent drug (CID 2733526), a different compound "
     "from OHT. Never pooled with OHT without a flag."),
    ("E2", "17beta-estradiol", "hormone", False, "5757", "The natural estrogen."),
    ("HRG_B1", "heregulin-beta1 (176-246)", "growth_factor", False, "",
     "A protein growth factor, so no PubChem CID. 10 nM per the S2 paper."),
    ("FULVESTRANT", "fulvestrant", "SERD", False, "104741", "Estrogen receptor degrader."),
    ("GDC0810", "GDC-0810 (brilanestrant)", "SERD", False, "56941241", "Experimental receptor degrader, S5."),
    ("GDC0927", "GDC-0927", "SERD", False, "87055263", "Experimental receptor degrader, S5."),
    ("GNE274", "GNE-274", "other", False, "162641011",
     "Non-degrading analog of GDC-0927 (S5 paper). Called 'G-03046274' in the S5 counts file. "
     "CID from PubChem name lookup."),
    ("ETHANOL", "ethanol", "vehicle", True, "702", "Solvent for vehicle controls in S1 (0.1%) and S3 (per papers)."),
    ("DMSO", "dimethyl sulfoxide", "vehicle", True, "679", "Solvent for vehicle controls in S5."),
]

# alias as written, compound_id, seen_in_study, seen_in_field
COMPOUND_ALIASES = [
    # 4-hydroxytamoxifen: the 9 spellings, plus the paper-only abbreviations
    ("TOT", "OHT", "S1", "title, characteristics, description"),
    ("trans-hydroxytamoxifen", "OHT", "S1", "extract protocol"),
    ("hydroxytamoxifen", "OHT", "S1", "series overall design"),
    ("Tam", "OHT", "S1", "paper only (paper_1.pdf abbreviation for trans-hydroxytamoxifen)"),
    ("4-hydroxytamoxifen", "OHT", "S3", "title, characteristics, source name"),
    ("OH-Tam", "OHT", "S3", "treatment protocol"),
    ("4-OH-tamoxifen", "OHT", "S4", "description, summary"),
    ("4-OH tamoxifen", "OHT", "S5", "characteristics (ligand)"),
    ("4OH-Tamoxifen", "OHT", "S5", "counts file header row 2"),
    ("4-OHT", "OHT", "S5", "paper only (paper_5.pdf)"),
    # tamoxifen of unstated form
    ("tamoxifen", "TAM_UNSPECIFIED", "S2", "characteristics (treatment)"),
    ("Tamoxifen", "TAM_UNSPECIFIED", "S2", "title, source name, description"),
    ("Tam", "TAM_UNSPECIFIED", "S2", "lab codes, e.g. TamR.E2_Tam.24h"),
    # estradiol
    ("E2", "E2", "S1", "title, characteristics"),
    ("17beta-estradiol", "E2", "S1", "series overall design, extract protocol"),
    ("E2", "E2", "S2", "title, characteristics, codes"),
    ("E2", "E2", "S5", "characteristics (ligand), counts file header"),
    # heregulin
    ("HRG", "HRG_B1", "S2", "title, characteristics, codes"),
    # receptor degraders
    ("Fulvestrant", "FULVESTRANT", "S5", "characteristics (ligand), counts file header"),
    ("GDC-0810", "GDC0810", "S5", "characteristics (ligand), counts file header"),
    ("GDC-0927", "GDC0927", "S5", "characteristics (ligand), counts file header"),
    ("GNE-274", "GNE274", "S5", "characteristics (ligand)"),
    ("G-03046274", "GNE274", "S5", "counts file header row 2"),
    # vehicles
    ("vehicle", "ETHANOL", "S1", "characteristics (solvent is ethanol per paper_1.pdf)"),
    ("veh", "ETHANOL", "S1", "title, characteristics"),
    ("vehicle", "ETHANOL", "S3", "title, source name (solvent is ethanol per paper_3.pdf)"),
    ("control (vehicle)", "ETHANOL", "S3", "characteristics (agent)"),
    ("DMSO", "DMSO", "S5", "characteristics (ligand)"),
    ("DMSO control", "DMSO", "S5", "counts file header row 2"),
]

# cell_line_id, preferred_name, is_derived, parent_cell_line_id, derived_by_study, derivation, disease, source_bank, notes
CELL_LINES = [
    ("CVCL_0031", "MCF-7", False, "", "", "", "Invasive breast carcinoma of no special type", "ATCC", "Used by all five studies."),
    ("CVCL_0553", "T-47D", False, "", "", "", "Invasive breast carcinoma of no special type", "ATCC", "S4, S5."),
    ("CVCL_0179", "BT-474", False, "", "", "", "Invasive breast carcinoma of no special type", "ATCC", "S4, S5."),
    ("CVCL_0588", "ZR-75-1", False, "", "", "", "Invasive breast carcinoma of no special type", "ATCC", "S4."),
    ("CVCL_1254", "HCC1500", False, "", "", "", "Breast ductal carcinoma", "ATCC", "S5."),
    ("CVCL_0619", "MDA-MB-330", False, "", "", "", "Invasive breast lobular carcinoma", "ATCC", "S5."),
    ("CVCL_1115", "CAMA-1", False, "", "", "", "Breast adenocarcinoma", "ATCC", "S5."),
    ("CVCL_0253", "EFM-19", False, "", "", "", "Breast ductal carcinoma", "DSMZ", "S5."),
    # S2 (RIKEN): six resistant clones, one chosen for the time course. The paper doesn't say which.
    ("S2_TAMR_SELECTED", "MCF-7 TamR, selected clone (RIKEN)", True, "CVCL_0031", "S2",
     "MCF-7 in 1 uM tamoxifen for 1 month, then 2 more months; one of six clones chosen", "", "lab-made",
     "Used for the time course. Identical to one of S2_TAMR_C1..C6, but the paper never says which."),
] + [
    (f"S2_TAMR_C{i}", f"MCF-7 TamR clone #{i} (RIKEN)", True, "CVCL_0031", "S2",
     "MCF-7 in 1 uM tamoxifen for 1 month, then 2 more months", "", "lab-made",
     "One of six resistant clones. Not a replicate of the others.")
    for i in range(1, 7)
] + [
    # S4 (Helsinki): resistant derivatives, 1 uM 4-OHT for 8-12 months, kept in drug afterwards
    ("S4_MCF7_TAM1", "MCF-7 Tam1 (Helsinki)", True, "CVCL_0031", "S4", "1 uM 4-OHT for 8-12 months", "", "lab-made", "Kept in 1 uM 4-OHT (paper)."),
    ("S4_T47D_TAM1", "T-47D Tam1 (Helsinki)", True, "CVCL_0553", "S4", "1 uM 4-OHT for 8-12 months", "", "lab-made", "Kept in 1 uM 4-OHT (paper)."),
    ("S4_T47D_TAM2", "T-47D Tam2 (Helsinki)", True, "CVCL_0553", "S4", "1 uM 4-OHT for 8-12 months", "", "lab-made", "Kept in 1 uM 4-OHT (paper)."),
    ("S4_ZR751_TAM1", "ZR-75-1 Tam1 (Helsinki)", True, "CVCL_0588", "S4", "1 uM 4-OHT for 8-12 months", "", "lab-made", "Kept in 1 uM 4-OHT (paper)."),
    ("S4_ZR751_TAM2", "ZR-75-1 Tam2 (Helsinki)", True, "CVCL_0588", "S4", "1 uM 4-OHT for 8-12 months", "", "lab-made", "Kept in 1 uM 4-OHT (paper)."),
    ("S4_BT474_TAM1", "BT-474 Tam1 (Helsinki)", True, "CVCL_0179", "S4", "1 uM 4-OHT for 8-12 months", "", "lab-made", "Kept in 1 uM 4-OHT (paper)."),
    ("S4_BT474_TAM2", "BT-474 Tam2 (Helsinki)", True, "CVCL_0179", "S4", "1 uM 4-OHT for 8-12 months", "", "lab-made", "Kept in 1 uM 4-OHT (paper)."),
]

# alias as written, cell_line_id, seen_in_study
CELL_LINE_ALIASES = [
    ("MCF-7 breast cancer cells", "CVCL_0031", "S1"),
    ("MCF-7", "CVCL_0031", "S1"),
    ("MCF-7", "CVCL_0031", "S2"),
    ("WT", "CVCL_0031", "S2"),
    ("MCF-7 breast cancer cell line", "CVCL_0031", "S3"),
    ("MCF-7cells", "CVCL_0031", "S3"),
    ("MCF-7", "CVCL_0031", "S4"),
    ("MCF-7", "CVCL_0031", "S5"),
    ("MCF7", "CVCL_0031", "S5"),
    ("T-47D", "CVCL_0553", "S4"),
    ("T-47D", "CVCL_0553", "S5"),
    ("T47D", "CVCL_0553", "S5"),
    ("BT-474", "CVCL_0179", "S4"),
    ("BT474", "CVCL_0179", "S4"),
    ("BT-474", "CVCL_0179", "S5"),
    ("ZR-75-1", "CVCL_0588", "S4"),
    ("HCC1500", "CVCL_1254", "S5"),
    ("MDA-MB-330", "CVCL_0619", "S5"),
    ("CAMA-1", "CVCL_1115", "S5"),
    ("Cama-1", "CVCL_1115", "S5"),
    ("EFM-19", "CVCL_0253", "S5"),
    ("TamR", "S2_TAMR_SELECTED", "S2"),
] + [(f"TamR#{i}", f"S2_TAMR_C{i}", "S2") for i in range(1, 7)] + [
    ("MCF-7 Tam1", "S4_MCF7_TAM1", "S4"),
    ("T-47D Tam1", "S4_T47D_TAM1", "S4"),
    ("T-47D Tam2", "S4_T47D_TAM2", "S4"),
    ("ZR-75-1 Tam1", "S4_ZR751_TAM1", "S4"),
    ("ZR-75-1 Tam2", "S4_ZR751_TAM2", "S4"),
    ("BT-474 Tam1", "S4_BT474_TAM1", "S4"),
    ("BT-474 Tam2", "S4_BT474_TAM2", "S4"),
]


def write(name, header, rows):
    with open(HERE / name, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for r in rows:
            w.writerow(["true" if v is True else "false" if v is False else v for v in r])
    print(f"{name}: {len(rows)} rows")


if __name__ == "__main__":
    write("compound.csv", ["compound_id", "preferred_name", "compound_class", "is_vehicle", "pubchem_cid", "notes"], COMPOUNDS)
    write("compound_alias.csv", ["alias", "compound_id", "seen_in_study", "seen_in_field"], COMPOUND_ALIASES)
    write("cell_line.csv", ["cell_line_id", "preferred_name", "is_derived", "parent_cell_line_id", "derived_by_study",
                            "derivation", "disease", "source_bank", "notes"], CELL_LINES)
    write("cell_line_alias.csv", ["alias", "cell_line_id", "seen_in_study"], CELL_LINE_ALIASES)
