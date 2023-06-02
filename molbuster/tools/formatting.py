import pandas as pd


def _create_long_output(df: pd.DataFrame) -> str:
    # return df.T.to_string(index=True, header=False)
    df = df.T
    cols = df.columns
    output = ""
    for col in cols:
        output += "--> " + " ".join(col) + "\n"
        output += df[[col]].to_string(index=True, header=False)
    return output


def _create_short_results(df: pd.DataFrame) -> str:
    results = df.copy()
    results.columns = results.columns.to_flat_index()
    columns = results.columns
    passes = results[columns].sum(axis=1)
    results["Passed tests"] = passes.apply(lambda x: f"passes ({x} / {len(columns)})")
    results["Pass"] = results[columns].all(axis=1)
    results.index.name = None
    results = results[["Passed tests"]]
    return results.to_string(index=True, header=False) + "\n"
