"""Package reviewed code/docs/test evidence without transient build files."""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--source-archive", required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    targets = ["winding_pipeline", "tests", "tools", "docs", "examples", "legacy_reference", "validation", "synthetic_example"]
    paths = [root / name for name in ("README.md", "pyproject.toml", "LICENSE", "NOTICE")]
    for folder in targets:
        paths.extend(p for p in (root / folder).rglob("*") if p.is_file() and "__pycache__" not in p.parts)
    paths = sorted(paths)
    hashes = {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    hashes["uploaded_source_archive_sha256"] = hashlib.sha256(Path(args.source_archive).read_bytes()).hexdigest()
    dest = Path(args.output).resolve()
    dest.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(dest, "x", compression=zipfile.ZIP_DEFLATED) as z:
        for p in paths:
            z.write(p, Path("awn_pipeline") / p.relative_to(root))
        z.writestr("awn_pipeline/RELEASE_SHA256.json", json.dumps(hashes, indent=2) + "\n")
    with zipfile.ZipFile(dest) as z:
        bad = z.testzip()
        if bad:
            raise RuntimeError(f"Bad ZIP member: {bad}")
        print(f"Verified ZIP: {len(z.namelist())} entries, {dest.stat().st_size} bytes")
    print(dest)


if __name__ == "__main__":
    main()
