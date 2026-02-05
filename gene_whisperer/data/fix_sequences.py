import pandas as pd  # pyright: ignore[reportMissingImports]

def clean_file(path, out_path):
    print(f"Cleaning {path} → {out_path}")
    df = pd.read_csv(path, sep="\t")

    # Assumes the sequence column is named `"sequence"` or `"seq"`
    seq_col = "sequence" if "sequence" in df.columns else "seq"

    def clean(seq):
        if not isinstance(seq, str):
            return seq
        # Remove spaces and uppercase
        seq = seq.replace(" ", "").replace("\t", "").strip().upper()
        # Keep ONLY A/C/G/T characters
        return "".join([b for b in seq if b in "ACGT"])

    df[seq_col] = df[seq_col].apply(clean)

    df.to_csv(out_path, sep="\t", index=False)
    print("Done.\n")

clean_file("stage2_train.tsv", "stage2_train.cleaned.tsv")
clean_file("stage2_val.tsv", "stage2_val.cleaned.tsv")
print("All cleaned!")