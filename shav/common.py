"""Small shared I/O helpers. Invalid or incomplete inputs raise explicitly."""
from __future__ import annotations

import csv
import hashlib
import json
import subprocess
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
SEGMENTS = ("L", "M", "S")
Row = dict[str, str]


def read_tsv(path: Path) -> list[Row]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Missing TSV header: {path}")
        rows = list(reader)
    if any(None in row or any(v is None for v in row.values()) for row in rows):
        raise ValueError(f"Malformed TSV: {path}")
    return rows


def write_tsv(path: Path, rows: Sequence[Row], fields: Sequence[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = list(fields) if fields is not None else list(dict.fromkeys(k for row in rows for k in row))
    if not columns:
        raise ValueError("A TSV requires explicit columns or at least one row")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    name = ""
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            header = line[1:].split()
            if not header:
                raise ValueError(f"Empty FASTA identifier: {path}")
            name = header[0]
            if not name or name in sequences:
                raise ValueError(f"Empty or duplicate FASTA identifier: {name}")
            sequences[name] = ""
        elif line.strip():
            if not name:
                raise ValueError(f"Sequence before FASTA header: {path}")
            sequences[name] += line.strip().upper()
    if not sequences or any(not sequence for sequence in sequences.values()):
        raise ValueError(f"Empty FASTA sequence: {path}")
    return sequences


def write_fasta(path: Path, sequences: Iterable[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        for name, sequence in sequences:
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def read_json(path: Path) -> Any:
    # JSON is a boundary format; callers validate its domain shape before use.
    return json.loads(path.read_text())


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(command: Sequence[str], *, stdout: Path | None = None, log: Path | None = None) -> None:
    if log:
        log.parent.mkdir(parents=True, exist_ok=True)
    if stdout:
        stdout.parent.mkdir(parents=True, exist_ok=True)
    with (log.open("w") if log else open("/dev/null", "w")) as err:
        if stdout:
            with stdout.open("w") as out:
                subprocess.run(list(command), check=True, stdout=out, stderr=err)
        else:
            subprocess.run(list(command), check=True, stdout=err, stderr=err)

