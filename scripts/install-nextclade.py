"""Install the pinned official Nextclade release for this platform."""
from __future__ import annotations

import argparse
import json
import platform
import subprocess
import urllib.request
from pathlib import Path

VERSION = "3.25.1"


def version(binary: Path) -> str:
    return subprocess.run([str(binary), "--version"], check=True, capture_output=True, text=True).stdout.strip()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true", help="Re-download and verify the official release")
    args = parser.parse_args()
    destination = Path(__file__).resolve().parents[1] / ".runtime/bin/nextclade"
    if destination.exists() and not args.force and version(destination) == f"nextclade {VERSION}":
        print(f"nextclade {VERSION} already installed")
        return
    suffixes = {("Darwin", "arm64"): "aarch64-apple-darwin", ("Darwin", "x86_64"): "x86_64-apple-darwin",
                ("Linux", "aarch64"): "aarch64-unknown-linux-gnu", ("Linux", "x86_64"): "x86_64-unknown-linux-gnu"}
    suffix = suffixes[(platform.system(), platform.machine())]
    request = urllib.request.Request(f"https://api.github.com/repos/nextstrain/nextclade/releases/tags/{VERSION}", headers={"User-Agent": "shav-phylogenomics"})
    with urllib.request.urlopen(request, timeout=30) as response:
        release = json.load(response)
    asset = next(a for a in release["assets"] if a["name"] == "nextclade-" + suffix)
    url = str(asset["browser_download_url"])
    if not url.startswith(f"https://github.com/nextstrain/nextclade/releases/download/{VERSION}/"):
        raise ValueError("Unexpected official release asset destination")
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(".download")
    # curl handles transient CDN failures and slow connections. Keep the installed
    # binary intact until the complete candidate passes size and version checks.
    subprocess.run(["curl", "--fail", "--location", "--silent", "--show-error", "--connect-timeout", "20",
                    "--max-time", "300", "--retry", "2", "--output", str(temporary), url], check=True)
    if temporary.stat().st_size != int(asset["size"]):
        raise ValueError("Incomplete Nextclade release download")
    temporary.chmod(0o755)
    actual = version(temporary)
    if actual != f"nextclade {VERSION}":
        raise ValueError(f"Wrong Nextclade version: {actual}")
    temporary.replace(destination)
    print(actual)


if __name__ == "__main__":
    main()
