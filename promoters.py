import gffutils
import pandas as pd

# Replace this with the path to your gene annotation GTF file.
gff_file = "homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20240230.gff"

db = gffutils.create_db(
    gff_file,
    dbfn="regulatory_features.db",
    force=True,
    keep_order=True,
    merge_strategy="merge",
    sort_attribute_values=True
)

promoter_records = []

for feature in db.all_features():
    if "promoter" in feature.featuretype.lower():
        promoter_records.append({
            "chrom": feature.chrom,
            "start": feature.start,
            "end": feature.end,
            "strand": feature.strand,
            "feature_type": feature.featuretype,
            "desc": feature.attributes['ID']
        })

df_promoters = pd.DataFrame(promoter_records)
