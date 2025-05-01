"""Functions for formatting the output of the PoseBusters tools."""

import pandas as pd


def _value_map(x):
    if isinstance(x, bool):
        return ".   " if x else "Fail"
    return x


def create_long_output(df: pd.DataFrame) -> str:
    """Create a detailed output string for the long PoseBusters reporting format.

    Args:
        df: Report dataframe.

    Returns:
        String with detailed report.
    """
    df = df.T
    df[df.columns[-1]] = df[df.columns[-1]].map(_value_map, na_action="ignore")
    cols = df.columns
    output = ""
    segments = []
    for col in cols:
        segments.append(
            "Long summary for " + " ".join(map(str, col)) + "\n" + df[[col]].to_string(index=True, header=False)
        )
    output += "\n\n".join(segments)
    return output + "\n"


def create_short_output(df: pd.DataFrame) -> str:
    """Create a one line output string for the short PoseBusters reporting format.

    Args:
        df: Report dataframe.

    Returns:
        String with one line report.
    """
    results = df.copy()
    results.columns = results.columns.to_flat_index()
    columns = results.columns
    passes = results[columns].sum(axis=1)
    results["Passed tests"] = passes.apply(lambda x: f"passes ({x} / {len(columns)})")
    results["Pass"] = results[columns].all(axis=1)
    results.index.name = None
    results = results[["Passed tests"]]
    return results.to_string(index=True, header=False) + "\n"
