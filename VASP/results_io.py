"""
results_io.py

Shared helper for reading/writing the pipeline's intermediate results file
(DPT_results.csv, and can be reused for other channels' result files too).

Data model: long format, one row per quantity:
    quantity | value | unit | source

This avoids the overwrite problem you get with a single wide row when
different notebooks (hydro, uniaxial, shear, elastic compliances, ...)
each contribute different quantities to the same file.

Usage
-----
from results_io import save_result, save_results, load_result, load_results

# writing (e.g. in DFT.ipynb, after computing a, Xiu, d, dEg_dP)
save_results(
    {"a": a, "Xi_u": Xiu, "|d|": abs(d), "dEg/dP": dEg_dP},
    unit={"a": "eV", "Xi_u": "eV", "|d|": "eV", "dEg/dP": "eV/GPa"},
    source="DFT.ipynb",
)

# reading (e.g. in hydro.ipynb)
a = load_result("a")

# reading everything as a dataframe, to eyeball the full table
df = load_results()
"""

import os
import ast
import json
import numpy as np
import pandas as pd

# Anchor to the folder this file lives in (e.g. VASP/), NOT the caller's cwd.
# This makes the default path work the same regardless of which notebook,
# or which subfolder depth, imports this module.
_MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_PATH = os.path.join(_MODULE_DIR, "DPT_results.csv")


def _to_storable(v):
    """
    Coerce an incoming value into something safe to round-trip through a
    single CSV cell:
      - plain scalar (int/float/np scalar)         -> python float
      - length-1 list/tuple/array/Series            -> unwrapped python float
      - longer list/tuple/array/Series               -> JSON string, e.g. "[1.2, 3.4]"
      - numeric string ("1.83")                      -> python float
      - list-looking string ("[1.2, 3.4]")           -> re-validated JSON string
    Raises ValueError (with the offending value shown) for anything else,
    rather than silently writing garbage that poisons the column later.
    """
    if isinstance(v, (pd.Series, pd.Index)):
        v = v.to_numpy()

    if isinstance(v, (list, tuple, np.ndarray)):
        arr = np.asarray(v, dtype=float)
        if arr.size == 1:
            return float(arr.reshape(-1)[0])
        # nested list -> shape is preserved through the round trip
        return json.dumps(arr.tolist())

    if isinstance(v, str):
        s = v.strip()
        try:
            return float(s)
        except ValueError:
            pass
        try:
            parsed = ast.literal_eval(s)
        except (ValueError, SyntaxError):
            raise ValueError(
                f"Could not interpret string value {v!r} as a number or list of numbers."
            )
        return _to_storable(parsed)

    try:
        return float(v)
    except (TypeError, ValueError):
        raise ValueError(
            f"Value {v!r} of type {type(v)} could not be converted to a number "
            f"or list of numbers."
        )


def _from_storable(raw):
    """Inverse of _to_storable: turn a stored cell back into a float or list of floats."""
    if isinstance(raw, (int, float, np.floating, np.integer)) and not isinstance(raw, bool):
        return float(raw)
    if isinstance(raw, str):
        s = raw.strip()
        if s.startswith("["):
            return json.loads(s)
        return float(s)
    raise ValueError(f"Could not interpret stored value {raw!r}.")


def _load_table(path=DEFAULT_PATH):
    """Load the results CSV, or an empty properly-shaped table if it doesn't exist yet."""
    if os.path.exists(path):
        return pd.read_csv(path)
    return pd.DataFrame(columns=["quantity", "value", "unit", "source"])


def save_result(quantity, value, unit=None, source=None, path=DEFAULT_PATH):
    """Upsert a single quantity into the results file."""
    save_results({quantity: value}, unit={quantity: unit}, source=source, path=path)


def save_results(values: dict, unit: dict = None, source: str = None, path=DEFAULT_PATH):
    """
    Upsert multiple quantities at once.

    values : dict of {quantity_name: value}
    unit   : optional dict of {quantity_name: unit_string}
    source : optional string (e.g. "DFT.ipynb") applied to all rows in this call
    """
    unit = unit or {}
    df = _load_table(path)

    new_rows = pd.DataFrame([
        {
            "quantity": q,
            "value": _to_storable(v),
            "unit": unit.get(q),
            "source": source,
        }
        for q, v in values.items()
    ])

    # drop any existing rows for these quantities, then append the fresh ones
    df = df[~df["quantity"].isin(new_rows["quantity"])]
    df = pd.concat([df, new_rows], ignore_index=True)

    df.to_csv(path, index=False)
    return df


def load_result(quantity, path=DEFAULT_PATH):
    """Return the value for a single quantity (raises KeyError if not found)."""
    df = _load_table(path)
    match = df.loc[df["quantity"] == quantity, "value"]
    if match.empty:
        available = list(df["quantity"])
        raise KeyError(
            f"'{quantity}' not found in {path}. Available quantities: {available}"
        )
    return _from_storable(match.iloc[0])


def load_array(quantity, path=DEFAULT_PATH):
    """Like load_result, but always returns a numpy array (shape preserved)."""
    return np.asarray(load_result(quantity, path=path), dtype=float)


def save_table(df, path):
    """
    Save a tidy table (one row per sample point, e.g. per strain step).

    Use this for SWEEP data -- arrays that share a common index -- rather than
    stuffing each array into a single cell of the long-format results file.
    """
    df = pd.DataFrame(df)
    df.to_csv(path, index=False)
    return df


def load_table(path):
    """Load a tidy table written by save_table."""
    if not os.path.exists(path):
        raise FileNotFoundError(f"No table at {path}")
    return pd.read_csv(path)


def load_results(quantities=None, path=DEFAULT_PATH, parse=True):
    """
    Return results as a dataframe.
    If `quantities` is given (a list), returns just those rows.
    If `parse` is True (default), the 'value' column is parsed per-cell into
    python floats/lists, so a single list-valued row can't make previously
    clean float rows look like strings. Set parse=False to see the raw
    on-disk cell contents (useful when debugging a corrupted file).
    """
    df = _load_table(path)
    if quantities is not None:
        df = df[df["quantity"].isin(quantities)].reset_index(drop=True)

    if parse and not df.empty:
        df = df.copy()
        df["value"] = df["value"].map(_from_storable)

    return df