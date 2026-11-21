#!/usr/bin/env python3
"""Convert the archived compact Group-B text instances to validated DDT files."""

import argparse
from dataclasses import dataclass
import math
from pathlib import Path
import re
import tempfile


ROOT = Path(__file__).resolve().parents[1]
BOUND = re.compile(
    r"^\s*\d+:\s*\[([-+0-9.eE]+)\s+([-+0-9.eE]+)\]\s*"
    r"\[([-+0-9.eE]+)\s+([-+0-9.eE]+)\]\s*$"
)
ROW = re.compile(r"^\s*\d+:\s*(.*?)\s*<\s*([-+0-9.eE]+)\s*$")
TERM = re.compile(r"([+-]?)([PV])(\d+)")


@dataclass(frozen=True)
class SparseRow:
    base: str
    terms: tuple[tuple[int, int], ...]
    rhs: str


@dataclass(frozen=True)
class Instance:
    first: int
    second: int
    bounds: tuple[tuple[str, str, str, str], ...]
    probability_rows: tuple[SparseRow, ...]
    value_rows: tuple[SparseRow, ...]


def parse_instance(path: Path) -> Instance:
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 2:
        raise ValueError(f"{path}: truncated header")
    try:
        first, second = map(int, lines[0].split())
        probability_count, value_count = map(int, lines[1].split())
    except ValueError as error:
        raise ValueError(f"{path}: invalid two-line header") from error
    variables = first + second
    expected = 2 + variables + probability_count + value_count
    if (first < 1 or second < 1 or probability_count < 0 or value_count < 0 or
            variables > 300 or len(lines) != expected):
        raise ValueError(
            f"{path}: dimensions/line count invalid (got {len(lines)}, expected {expected})"
        )

    bounds = []
    for number, text in enumerate(lines[2:2 + variables], start=3):
        match = BOUND.fullmatch(text)
        if not match:
            raise ValueError(f"{path}:{number}: invalid bound row")
        values = match.groups()
        numeric = tuple(map(float, values))
        if not (0 <= numeric[0] <= numeric[1] <= 1 and
                0 <= numeric[2] <= numeric[3] <= 1):
            raise ValueError(f"{path}:{number}: bounds must lie in [0,1]")
        bounds.append(values)

    def parse_rows(texts: list[str], expected_base: str, start: int) -> tuple[SparseRow, ...]:
        result = []
        for offset, text in enumerate(texts):
            match = ROW.fullmatch(text)
            if not match:
                raise ValueError(f"{path}:{start + offset}: invalid linear row")
            expression, rhs = match.groups()
            if not math.isfinite(float(rhs)):
                raise ValueError(f"{path}:{start + offset}: non-finite rhs")
            matches = list(TERM.finditer(expression))
            if not matches or "".join(item.group(0) for item in matches) != expression:
                raise ValueError(f"{path}:{start + offset}: invalid expression")
            terms = []
            used = set()
            for item in matches:
                sign, base, raw_index = item.groups()
                index = int(raw_index)
                if base != expected_base or index < 1 or index > variables or index in used:
                    raise ValueError(f"{path}:{start + offset}: invalid variable in row")
                used.add(index)
                terms.append((index, -1 if sign == "-" else 1))
            result.append(SparseRow(expected_base, tuple(terms), rhs))
        return tuple(result)

    row_start = 2 + variables
    p_text = lines[row_start:row_start + probability_count]
    v_text = lines[row_start + probability_count:]
    return Instance(
        first, second, tuple(bounds),
        parse_rows(p_text, "P", row_start + 1),
        parse_rows(v_text, "V", row_start + probability_count + 1),
    )


def address(index: int, first: int) -> tuple[int, int]:
    return (1, index) if index <= first else (2, index - first)


def encode(path: Path, instance: Instance) -> str:
    lines = [
        "8.8.4",
        path.stem,
        "1",
        "0",
        f"2 {instance.first} {instance.second}",
        "alternative_A",
        "alternative_B",
        str(len(instance.bounds) + len(instance.probability_rows)),
    ]
    for index, bound in enumerate(instance.bounds, start=1):
        alternative, consequence = address(index, instance.first)
        lines.append(f"1 {alternative} {consequence} 1 {bound[0]} {bound[1]}")
    for row in instance.probability_rows:
        fields = ["3", str(len(row.terms))]
        for index, coefficient in row.terms:
            alternative, consequence = address(index, instance.first)
            fields.extend((str(alternative), str(consequence), str(coefficient)))
        fields.append(row.rhs)
        lines.append(" ".join(fields))

    lines.extend(("0 1", str(len(instance.bounds) + len(instance.value_rows))))
    for index, bound in enumerate(instance.bounds, start=1):
        alternative, consequence = address(index, instance.first)
        lines.append(f"1 {alternative} {consequence} 1 {bound[2]} {bound[3]}")
    for row in instance.value_rows:
        fields = ["3", str(len(row.terms))]
        for index, coefficient in row.terms:
            alternative, consequence = address(index, instance.first)
            fields.extend((str(alternative), str(consequence), str(coefficient)))
        fields.append(row.rhs)
        lines.append(" ".join(fields))
    lines.extend(("0", "0"))
    return "\n".join(lines) + "\n"


def atomic_write(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent,
                                     prefix=path.name + ".", delete=False) as handle:
        temporary = Path(handle.name)
        handle.write(content)
    temporary.chmod(0o644)
    temporary.replace(path)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path, nargs="?", default=ROOT / "data/B")
    parser.add_argument("--output-root", type=Path,
                        help="mirror files below this directory; default writes beside each .m")
    parser.add_argument("--check", action="store_true",
                        help="validate sources and existing outputs without writing")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    source = args.source.resolve()
    source_base = source if source.is_dir() else source.parent
    files = sorted(source.rglob("*.m")) if source.is_dir() else [source]
    if not files:
        parser.error(f"no .m instances found below {source}")

    written = checked = 0
    for input_path in files:
        instance = parse_instance(input_path)
        content = encode(input_path, instance)
        output = ((args.output_root.resolve() / input_path.relative_to(source_base)).with_suffix(".ddt")
                  if args.output_root else input_path.with_suffix(".ddt"))
        if args.check:
            if not output.exists() or output.read_text(encoding="utf-8") != content:
                raise ValueError(f"missing or stale output: {output}")
            checked += 1
        else:
            if output.exists() and not args.overwrite:
                if output.read_text(encoding="utf-8") != content:
                    raise ValueError(f"refusing to replace differing output: {output}")
                checked += 1
            else:
                atomic_write(output, content)
                written += 1
    print(f"converted={written} already_current={checked} total={len(files)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
