"""Install the pinned official Nextclade release for this platform."""
from __future__ import annotations

import json
import platform
import subprocess
import urllib.request
from pathlib import Path

VERSION = "3.25.1"
root = Path(__file__).resolve().parents[1]
suffixes = {("Darwin", "arm64"): "aarch64-apple-darwin", ("Darwin", "x86_64"): "x86_64-apple-darwin",
            ("Linux", "aarch64"): "aarch64-unknown-linux-gnu", ("Linux", "x86_64"): "x86_64-unknown-linux-gnu"}
suffix = suffixes[(platform.system(), platform.machine())]
asset = "nextclade-" + suffix
request = urllib.request.Request(f"https://api.github.com/repos/nextstrain/nextclade/releases/tags/{VERSION}", headers={"User-Agent": "shav-phylogenomics"})
with urllib.request.urlopen(request, timeout=30) as response:
    release = json.load(response)
url = next(a["browser_download_url"] for a in release["assets"] if a["name"] == asset)
if not url.startswith(f"https://github.com/nextstrain/nextclade/releases/download/{VERSION}/"):
    raise ValueError("Unexpected official release asset destination")
destination = root / ".runtime/bin/nextclade"
destination.parent.mkdir(parents=True, exist_ok=True)
temporary = destination.with_suffix(".download")
with urllib.request.urlopen(url, timeout=120) as response:
    temporary.write_bytes(response.read())
temporary.chmod(0o755)
actual = subprocess.run([str(temporary), "--version"], check=True, capture_output=True, text=True).stdout.strip()
if actual != f"nextclade {VERSION}":
    raise ValueError(f"Wrong Nextclade version: {actual}")
temporary.replace(destination)
print(actual)
