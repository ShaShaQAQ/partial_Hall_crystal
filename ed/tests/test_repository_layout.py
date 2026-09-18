from pathlib import Path
import subprocess


ROOT = Path(__file__).resolve().parents[2]


def test_long_lived_modules_have_explicit_top_level_owners():
    expected = [
        ROOT / "ed" / "shared",
        ROOT / "ed" / "cases",
        ROOT / "ed" / "projected",
        ROOT / "dmrg",
        ROOT / "reports" / "fqahc_results",
        ROOT / "results",
    ]

    assert all(path.is_dir() for path in expected)
    assert not (ROOT / "shared").exists()
    assert not (ROOT / "projected_ed").exists()
    assert not any(ROOT.glob("case[0-9]*"))
    assert not (ROOT / "dmrg" / "report").exists()


def test_ed_sources_do_not_reference_the_legacy_root_shared_directory():
    offenders = []
    for path in (ROOT / "ed").rglob("*.jl"):
        if 'include("../shared/' in path.read_text():
            offenders.append(path.relative_to(ROOT).as_posix())

    assert offenders == []


def test_large_runtime_artifacts_are_ignored():
    probes = [
        "results/example/wavefunction.jld2",
        "ed/cases/example/output/run.out",
        "ed/cases/example/output/run.err",
    ]
    process = subprocess.run(
        ["git", "check-ignore", "--stdin"],
        cwd=ROOT,
        input="\n".join(probes) + "\n",
        text=True,
        capture_output=True,
        check=False,
    )

    assert process.returncode == 0
    ignored = set(process.stdout.splitlines())
    assert ignored == set(probes)
