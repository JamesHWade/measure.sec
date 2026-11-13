"""Independent fixed-fixture SEC unit/alignment oracle, Python standard library only."""

import bisect
import hashlib
import json
import math
import os
import platform
import re
import stat
import sys
from pathlib import Path

if not __debug__ or os.environ.get("PYTHONOPTIMIZE"):
    raise SystemExit(
        "Run without Python optimization; verification assertions must stay enabled."
    )

PROFILE = "sec-elution-operations/1"
FROZEN = {
    "fixtures.json": "dbd3c3e15015f81264a83743b5989da1fe0e43712fda71e3a823ab46e98b6680",
    "declaration.json": "26f2810af74c8a0e80eda4b0baaf29192f3d12e6c6adccf2d6eb68335a7e1316",
}
COUNTS = {"curves": 0, "locations": 0, "signal_values": 0, "diagnostic_rows": 0}


def unique_object(pairs):
    answer = {}
    for key, value in pairs:
        if key in answer:
            raise ValueError("Duplicate JSON field: " + key)
        answer[key] = value
    return answer


def invalid_constant(value):
    raise ValueError("Non-JSON numeric token: " + value)


def finite_float(value):
    result = float(value)
    if not math.isfinite(result):
        raise ValueError("JSON number exceeds finite binary64 range: " + value)
    return result


def read_json(path, max_bytes=32 * 1024 * 1024):
    info = path.lstat()
    assert stat.S_ISREG(info.st_mode) and info.st_size <= max_bytes, path
    with path.open("rb") as stream:
        raw = stream.read(max_bytes + 1)
    assert len(raw) <= max_bytes, path
    return json.loads(
        raw.decode("utf-8"),
        object_pairs_hook=unique_object,
        parse_constant=invalid_constant,
        parse_float=finite_float,
    ), hashlib.sha256(raw).hexdigest()


def number(value):
    assert isinstance(value, (int, float)) and not isinstance(value, bool)
    assert math.isfinite(value)
    return float(value)


def close(actual, expected):
    assert math.isclose(number(actual), expected, rel_tol=1e-12, abs_tol=1e-12), (
        actual,
        expected,
    )


def check_curve(observed, locations, values):
    assert isinstance(observed, dict) and set(observed) == {"location", "value"}
    assert isinstance(observed["location"], list) and isinstance(
        observed["value"], list
    )
    assert (
        len(observed["location"])
        == len(locations)
        == len(observed["value"])
        == len(values)
    )
    for actual, expected in zip(observed["location"], locations):
        assert number(actual) == expected, (actual, expected)
    for actual, expected in zip(observed["value"], values):
        assert isinstance(actual, dict) and set(actual) == {"state", "number"}
        if expected is None:
            assert actual == {"state": "na", "number": None}, actual
        else:
            assert actual["state"] == "number"
            close(actual["number"], expected)
    COUNTS["curves"] += 1
    COUNTS["locations"] += len(locations)
    COUNTS["signal_values"] += len(values)


def raw_area(locations, values):
    return math.fsum(
        abs(right - left) * (a + b) / 2
        for left, right, a, b in zip(locations, locations[1:], values, values[1:])
    )


def volume_factor(unit, flow):
    if unit == "min":
        return flow
    if unit == "s":
        return flow / 60
    return 1 if unit == "mL" else 0.001


def mass_expected(fixture, case):
    x = [time * case["axis_per_minute"] for time in fixture["time_min"]]
    y = fixture["signal"][:]
    if case["reverse"]:
        x.reverse()
        y.reverse()
    positive = [max(value, 0) for value in y]
    negative = [max(-value, 0) for value in y]
    pos_area, neg_area = raw_area(x, positive), raw_area(x, negative)
    k = volume_factor(case["elution_units"], case["flow_ml_min"])
    denominator = pos_area * k
    scale = 0.001 if case["concentration_units"] == "g/mL" else 1
    result = [
        value * fixture["injection_mass_mg"] / denominator * scale for value in positive
    ]
    diagnostics = {
        "negative_points": sum(value < 0 for value in y),
        "points": len(y),
        "positive_area": pos_area,
        "negative_area": neg_area,
        "excluded_fraction": neg_area / (pos_area + neg_area),
    }
    recovered = raw_area(x, result) * k * (1000 if scale == 0.001 else 1)
    close(recovered, fixture["injection_mass_mg"])
    return x, result, diagnostics


def delta_in_axis(case):
    value = case["delay_value"]
    unit = case["elution_units"]
    flow = case["flow_ml_min"]
    if case["delay_kind"] == "time":
        if unit == "min":
            return value
        if unit == "s":
            return value * 60
        return value * flow * (1000 if unit == "uL" else 1)
    if unit == "mL":
        return value
    if unit == "uL":
        return value * 1000
    return value / flow * (60 if unit == "s" else 1)


def interpolate(locations, values, queries):
    pairs = sorted(zip(locations, values))
    xs = [pair[0] for pair in pairs]
    ys = [pair[1] for pair in pairs]
    answer = []
    for query in queries:
        i = bisect.bisect_left(xs, query)
        if i < len(xs) and xs[i] == query:
            answer.append(ys[i])
        elif i == 0 or i == len(xs):
            answer.append(None)
        elif ys[i - 1] is None or ys[i] is None:
            answer.append(None)
        else:
            weight = (query - xs[i - 1]) / (xs[i] - xs[i - 1])
            answer.append(ys[i - 1] + weight * (ys[i] - ys[i - 1]))
    return answer


def delay_expected(fixture, case):
    x = [time * case["axis_per_minute"] for time in fixture["time_min"]]
    y, reference = fixture["signal"][:], fixture["reference"][:]
    if case["reverse"]:
        x.reverse()
        y.reverse()
        reference.reverse()
    delta = delta_in_axis(case)
    if delta == 0:
        return x, y, reference
    if case["method"] == "interpolate":
        shifted = interpolate(x, y, [location + delta for location in x])
    else:
        index = round(delta / (x[1] - x[0]))
        shifted = [
            y[i + index] if 0 <= i + index < len(y) else None for i in range(len(y))
        ]
    return x, shifted, reference


def indexed_cases(observed, declared):
    assert isinstance(observed, list) and len(observed) == len(declared)
    specifications = {case["id"]: case for case in declared}
    assert len(specifications) == len(declared)
    seen = set()
    for case in observed:
        assert isinstance(case, dict)
        name = case["id"]
        assert name in specifications and name not in seen, name
        seen.add(name)
        yield case, specifications[name]
    assert seen == set(specifications)


def main():
    assert len(sys.argv) == 2, "Provide the qualification output root."
    output = Path(sys.argv[1])
    here = Path(__file__).resolve().parent
    inputs = here / "qualification-input"
    if not inputs.is_dir():
        inputs = here
    fixtures, fixture_sha = read_json(inputs / "fixtures.json")
    declaration, declaration_sha = read_json(inputs / "declaration.json")
    hashes = {"fixtures.json": fixture_sha, "declaration.json": declaration_sha}
    assert hashes == FROZEN and fixtures["profile"] == declaration["profile"] == PROFILE
    observed, result_sha = read_json(output / "results.json")
    runtime, runtime_sha = read_json(output / "runtime.json")
    status, status_sha = read_json(output / "run-status.json")
    assert observed["profile"] == PROFILE and observed["input_files_sha256"] == FROZEN
    assert (
        status["complete"] is True
        and status["stage"] == "complete"
        and status["error"] is None
    )
    assert isinstance(status["updated_utc"], str) and status["updated_utc"]
    assert isinstance(runtime["created_utc"], str) and runtime["created_utc"]
    session = runtime["session_info"]
    assert (
        isinstance(session, list)
        and session
        and all(isinstance(line, str) for line in session)
    )
    assert session[0] == observed["R"]
    versions = runtime["dependency_versions"]
    for dependency in ["measure", "measure.sec", "recipes", "jsonlite", "digest"]:
        assert isinstance(versions[dependency], str) and versions[dependency]
    packages = observed["installed_packages"]
    assert set(packages) == {"measure", "measure.sec"}
    for package, identity in packages.items():
        assert isinstance(identity["path"], str) and identity["path"]
        assert identity["namespace_path"] == identity["path"]
        assert identity["version"] == versions[package]
    head = observed["tool_checkout_head"]
    assert (
        head is None
        or isinstance(head, str)
        and re.fullmatch(r"(?:[0-9a-f]{40}|[0-9a-f]{64})", head)
    )

    mass_fixtures = {f["id"]: f for f in fixtures["mass_traces"]}
    for actual, case in indexed_cases(observed["mass"], declaration["mass"]):
        assert actual["concentration_unit"] == case["concentration_units"]
        members = actual["members"]
        assert isinstance(members, list) and len(members) == len(case["bake_order"])
        expected_diagnostics = []
        for row, (member, name) in enumerate(zip(members, case["bake_order"]), 1):
            assert member["fixture_id"] == name and number(member["row"]) == row
            fixture = mass_fixtures[name]
            close(member["injection_mass_mg"], fixture["injection_mass_mg"])
            x, y, diagnostics = mass_expected(fixture, case)
            check_curve(member["curve"], x, y)
            expected_diagnostics.append({"row": row, **diagnostics})
        condition = actual["condition"]
        assert (
            condition["column"] == "ri"
            and condition["elution_units"] == case["elution_units"]
        )
        diagnostics = condition["diagnostics"]
        assert isinstance(diagnostics, list) and len(diagnostics) == len(
            expected_diagnostics
        )
        for actual_row, expected in zip(diagnostics, expected_diagnostics):
            assert set(actual_row) == set(expected)
            for key in expected:
                if key in ["row", "negative_points", "points"]:
                    assert number(actual_row[key]) == expected[key]
                else:
                    close(actual_row[key], expected[key])
            COUNTS["diagnostic_rows"] += 1

    delay_fixtures = {f["id"]: f for f in fixtures["delay_traces"]}
    for actual, case in indexed_cases(observed["delay"], declaration["delay"]):
        x, y, reference = delay_expected(delay_fixtures[case["fixture"]], case)
        check_curve(actual["target"], x, y)
        check_curve(actual["reference"], x, reference)
        assert actual["concentration_unit"] == "mg/mL"

    for actual, case in indexed_cases(observed["noop"], declaration["noop"]):
        fixture = fixtures["noop_trace"]
        check_curve(actual["curve"], fixture["time_min"], fixture["signal"])
        assert actual["concentration_unit"] == fixture["concentration_unit"]

    for actual, case in indexed_cases(observed["lals"], declaration["lals"]):
        fixture = fixtures["lals"]
        x = [t * case["axis_per_minute"] for t in fixture["time_min"]]
        spec = {
            "elution_units": case["elution_units"],
            "concentration_units": case["concentration_units"],
            "flow_ml_min": fixture["flow_ml_min"],
            "axis_per_minute": case["axis_per_minute"],
            "reverse": False,
        }
        cfixture = {
            "time_min": fixture["time_min"],
            "signal": fixture["ri_signal"],
            "injection_mass_mg": fixture["injection_mass_mg"],
        }
        _, concentration, _ = mass_expected(cfixture, spec)
        queries = [
            location + fixture["delay_time_min"] * case["axis_per_minute"]
            for location in x
        ]
        corrected = interpolate(x, concentration, queries)
        c_g_ml = [
            value / 1000
            if case["concentration_units"] == "mg/mL" and value is not None
            else value
            for value in corrected
        ]
        mw = [
            signal * fixture["calibration_constant"] / (fixture["optical_constant"] * c)
            if c is not None
            else None
            for signal, c in zip(fixture["scattering_signal"], c_g_ml)
        ]
        for value in mw[:-1]:
            close(value, fixture["expected_molar_mass"])
        assert mw[-1] is None
        check_curve(actual["concentration"], x, corrected)
        check_curve(actual["scattering"], x, fixture["scattering_signal"])
        check_curve(actual["mw"], x, mw)
        assert actual["concentration_unit"] == case["concentration_units"]

    invalid = observed["invalid_cases"]
    assert isinstance(invalid, list) and len(invalid) == len(
        declaration["invalid_cases"]
    )
    assert {case["id"] for case in invalid} == set(declaration["invalid_cases"])
    for case in invalid:
        assert case["error"] is True
        assert isinstance(case["message"], str) and case["message"].strip()
    receipt = {
        "profile": PROFILE,
        "independent_reader": "Python standard-library sampled-area and alignment oracle/1",
        "python_version": platform.python_version(),
        "input_files_sha256": hashes,
        "results_sha256": result_sha,
        "runtime_sha256": runtime_sha,
        "run_status_sha256": status_sha,
        "counts": COUNTS,
        "mass_cases": len(declaration["mass"]),
        "delay_cases": len(declaration["delay"]),
        "noop_cases": len(declaration["noop"]),
        "lals_cases": len(declaration["lals"]),
        "expected_rejections": len(invalid),
        "all_declared_cases_match": True,
        "scope": declaration["scope"],
    }
    with (output / "independent-python-receipt.json").open(
        "x", encoding="utf-8"
    ) as stream:
        json.dump(receipt, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(
        json.dumps(
            {
                k: v
                for k, v in receipt.items()
                if k.endswith("cases")
                or k in ["counts", "expected_rejections", "all_declared_cases_match"]
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
