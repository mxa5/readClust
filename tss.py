import pandas as pd

GTF_PATH = "Homo_sapiens.GRCh38.113.gtf"

gtf = pd.read_csv(
    GTF_PATH,
    sep="\t",
    comment="#",
    header=None,
    names=["Chromosome", "Source", "Feature", "Start", "End", "Score", "Strand", "Frame", "Attributes"]
)

genes = gtf[gtf["Feature"] == "gene"].copy()

genes["gene_id"] = genes["Attributes"].str.extract(r'gene_id "([^"]+)"')
genes["gene_name"] = genes["Attributes"].str.extract(r'gene_name "([^"]+)"')

genes["gene_name"] = genes["gene_name"].fillna(genes["gene_id"])

genes["TSS"] = genes.apply(
    lambda row: row["Start"] if row["Strand"] == "+" else row["End"],
    axis=1
)

tss_df = genes[["Chromosome", "TSS", "gene_name", "Strand"]].copy()
tss_df.rename(columns={"TSS": "Start"}, inplace=True)
tss_df["End"] = tss_df["Start"] + 1

tss_df = tss_df[["Chromosome", "Start", "End", "gene_name", "Strand"]].reset_index(drop=True)
