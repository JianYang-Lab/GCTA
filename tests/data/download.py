#!/usr/bin/env python3

import argparse
import hashlib
import json
from pathlib import Path
import tempfile
import urllib.request

DATA_DIR = Path(__file__).resolve().parent
MANIFEST_PATH = DATA_DIR / "manifest.json"


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def source_url(manifest, path):
    return (
        f"{manifest['source_repository']}/raw/{manifest['source_commit']}/"
        f"{manifest['source_directory']}/{path.name}"
    )


def validate_plink(prefix):
    fam_path = prefix.with_suffix(".fam")
    bim_path = prefix.with_suffix(".bim")
    bed_path = prefix.with_suffix(".bed")
    with fam_path.open(encoding="utf-8") as stream:
        sample_count = sum(1 for line in stream if line.strip())
    with bim_path.open(encoding="utf-8") as stream:
        marker_count = sum(1 for line in stream if line.strip())
    with bed_path.open("rb") as stream:
        assert stream.read(3) == b"\x6c\x1b\x01", f"Invalid BED header: {bed_path}"
    expected_size = 3 + marker_count * ((sample_count + 3) // 4)
    actual_size = bed_path.stat().st_size
    assert actual_size == expected_size, (
        f"Invalid BED size: {bed_path}: {actual_size} != {expected_size}"
    )


def main():
    parser = argparse.ArgumentParser(description="Download and verify GCTA test data")
    parser.add_argument("--check", action="store_true", help="verify local files without downloading")
    args = parser.parse_args()

    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    for item in manifest["files"]:
        path = DATA_DIR / item["path"]
        if path.exists():
            actual = sha256(path)
            assert actual == item["sha256"], f"Checksum mismatch: {path}"
            continue
        assert not args.check, f"Missing test data: {path}"
        path.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as stream:
            temporary_path = Path(stream.name)
        try:
            with urllib.request.urlopen(source_url(manifest, path), timeout=60) as response:
                temporary_path.write_bytes(response.read())
            actual = sha256(temporary_path)
            assert actual == item["sha256"], f"Downloaded checksum mismatch: {path}"
            temporary_path.replace(path)
        finally:
            temporary_path.unlink(missing_ok=True)

    for prefix in (DATA_DIR / "gcta_legacy" / "test", DATA_DIR / "gcta_legacy" / "toy"):
        validate_plink(prefix)

    print(f"Verified {len(manifest['files'])} test data files")


if __name__ == "__main__":
    main()
