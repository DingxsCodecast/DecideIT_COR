"""Reconstruct pairwise Upper Tisza comparisons from FLAIRS 2003 Tables 1--4.

Source: M. Danielson, L. Ekenberg, J. Johansson (2003), "Decision
Evaluation of Three Flood Management Strategies", FLAIRS, pp. 491--495.
https://cdn.aaai.org/FLAIRS/2003/Flairs03-096.pdf

This is a transparent reanalysis of the published simulation intervals, not
the unpublished stakeholder-weighted model or a reproduction of its figure.
"""

from decimal import Decimal, getcontext
from pathlib import Path

getcontext().prec = 35
D = Decimal
EVENT_PROBABILITIES = [(D(a), D(b)) for a, b in [
    (".858", ".906"), (".034", ".052"), (".021", ".032"),
    (".028", ".041"), (".011", ".017"),
]]

# Every tuple is (minimum, maximum) net income in million HUF. The paper's
# printed tables generally display maximum before minimum; this canonical
# ordering has been checked against the actual PDF rather than OCR text.
NET_INCOME = {
    "government": {
        1: [("0", "0"), ("-11.2", "-7.50"), ("-147", "-97.8"),
            ("-273", "-182"), ("-953", "-636")],
        2: [("0", "0"), ("-5.62", "-3.75"), ("-73.3", "-48.9"),
            ("-136", "-90.9"), ("-477", "-318")],
        3: [("1.77", "2.66"), ("-8.59", "-5.73"), ("-144", "-96"),
            ("-270", "-180"), ("-951", "-634")],
    },
    "insurers": {
        1: [("1.82", "2.73"), ("-4.72", "-3.15"),
            ("-65.4", "-43.6"), ("-115", "-76.8"), ("-376", "-251")],
        2: [("1.98", "2.96"), ("-4.89", "-3.26"),
            ("-37.6", "-25.1"), ("-68.5", "-45.7"), ("-254", "-170")],
        3: [("0", "0")] * 5,
    },
    "property_owners": {
        1: [("-2.73", "-1.82"), ("-21.5", "-14.3"),
            ("-277", "-185"), ("-521", "-347"), ("-1850", "-1230")],
        2: [("-2.96", "-1.98"), ("-27", "-18"),
            ("-378", "-252"), ("-704", "-469"), ("-2450", "-1630")],
        3: [("-2.66", "-1.77"), ("-28.9", "-19.3"),
            ("-345", "-230"), ("-639", "-426"), ("-2230", "-1480")],
    },
}
ORIGINAL_MIN = D("-2450")
ORIGINAL_MAX = D("2.96")
SCALE = ORIGINAL_MAX - ORIGINAL_MIN


def canonical_intervals(perspective, scenario):
    return [(D(lo), D(hi)) for lo, hi in NET_INCOME[perspective][scenario]]


def fmt(number):
    # Fixed precision makes .ddt output stable across Python versions.
    return format(number, ".16f")


def verify_source_transcription():
    assert sum(lo for lo, _ in EVENT_PROBABILITIES) <= 1
    assert sum(hi for _, hi in EVENT_PROBABILITIES) >= 1
    all_values = []
    for perspective in NET_INCOME.values():
        assert set(perspective) == {1, 2, 3}
        for intervals in perspective.values():
            assert len(intervals) == 5
            for lo, hi in [(D(a), D(b)) for a, b in intervals]:
                assert lo <= hi
                all_values += [lo, hi]
    assert min(all_values) == ORIGINAL_MIN
    assert max(all_values) == ORIGINAL_MAX


def build_case(perspective, first, second):
    assert first < second and perspective in NET_INCOME
    probability_rows = []
    for alt in (1, 2):
        for k, (lo, hi) in enumerate(EVENT_PROBABILITIES, 1):
            probability_rows.append(f"1 {alt} {k} 1 {fmt(lo)} {fmt(hi)}")
    # The five flood events are common to the alternatives. Two general
    # linear inequalities encode each cross-alternative equality.
    for k in range(1, 6):
        probability_rows += [
            f"3 2 1 {k} 1 2 {k} -1 0",
            f"3 2 1 {k} -1 2 {k} 1 0",
        ]
    value_rows = []
    for alt, scenario in [(1, first), (2, second)]:
        for k, (lo, hi) in enumerate(canonical_intervals(perspective, scenario), 1):
            lo01 = (lo - ORIGINAL_MIN) / SCALE
            hi01 = (hi - ORIGINAL_MIN) / SCALE
            assert D(0) <= lo01 <= hi01 <= D(1)
            value_rows.append(f"1 {alt} {k} 1 {fmt(lo01)} {fmt(hi01)}")
    lines = [
        "8.7.8", "Upper Tisza flood insurance; FLAIRS 2003 Tables 1-4",
        "5 common flood events; perspective=" + perspective,
        "Net incomes affine-normalized globally to [0,1]",
        "2 5 5", f"Scenario_{first}", f"Scenario_{second}",
        str(len(probability_rows)), *probability_rows,
        "values", str(len(value_rows)), *value_rows,
    ]
    return "\n".join(lines) + "\n"


def main():
    verify_source_transcription()
    out = Path(__file__).resolve().parent / "cases"
    out.mkdir(exist_ok=True)
    for perspective in NET_INCOME:
        for first, second in ((1, 2), (1, 3), (2, 3)):
            name = f"upper_tisza_{perspective}_s{first}_s{second}.ddt"
            (out / name).write_text(build_case(perspective, first, second),
                                    encoding="utf-8")
            print(out / name)


if __name__ == "__main__":
    main()
