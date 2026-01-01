"""Recover the audited CC0 PAIR-MaP subset and versioned reference sequence."""
from pathlib import Path
import hashlib
import io
import json
import urllib.request
import zipfile

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data/pairmap"


def sha256(content):
    return hashlib.sha256(content).hexdigest()


def main():
    manifest = json.loads((DATA / "source_manifest.json").read_text())
    payload = urllib.request.urlopen(manifest["source_url"], timeout=120).read()
    if sha256(payload) != manifest["archive_sha256"]:
        raise ValueError("PAIR-MaP archive checksum changed")
    archive = zipfile.ZipFile(io.BytesIO(payload))
    for record in manifest["files"]:
        relative = Path(record["path"])
        member = relative.relative_to("data/pairmap/source").as_posix()
        content = archive.read(member)
        if sha256(content) != record["sha256"]:
            raise ValueError(f"Source checksum mismatch: {member}")
        destination = ROOT / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(content)
    reference = DATA / "U00096.2.gb"
    if reference.exists() and sha256(reference.read_bytes()) == manifest["reference_genbank_sha256"]:
        return
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
           "db=nuccore&id=U00096.2&rettype=gbwithparts&retmode=text")
    content = urllib.request.urlopen(url, timeout=120).read()
    if sha256(content) != manifest["reference_genbank_sha256"]:
        raise ValueError("Reference annotation response differs from frozen GenBank bytes")
    reference.write_bytes(content)


if __name__ == "__main__":
    main()
