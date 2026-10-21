"""Project release gate; assigning the author's license is an explicit decision."""
import os
import re
from pathlib import Path

root = Path(__file__).resolve().parents[1]
metadata = (root/"pyproject.toml").read_text(encoding="utf-8")
version = re.search(r'^version\s*=\s*"([^"]+)"', metadata, re.M).group(1)
if os.environ.get("GITHUB_REF_TYPE") != "tag" or os.environ.get("GITHUB_REF_NAME") != "v"+version:
    raise SystemExit("Dispatch this workflow from the matching v"+version+" release tag.")
if not (root/"LICENSE").is_file() or not re.search(r'^license\s*=', metadata, re.M):
    raise SystemExit("Owner must select LICENSE and pyproject license metadata before release.")
print("Release tag, version and license metadata are present.")
