"""Build one portable HTML manual with local navigation and pre-rendered math.

Requires the docs+symbolic extras and Node.js at build time. The delivered
HTML needs neither JavaScript nor companion files and can be opened offline.
"""
from __future__ import annotations

import argparse
import base64
import json
import mimetypes
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from urllib.parse import unquote, urlsplit

from bs4 import BeautifulSoup
import sphinxcontrib.katex


ROOT = Path(__file__).resolve().parents[1]


def asset(root, reference):
    """Resolve an existing generated asset within this documentation build."""
    path = (root/unquote(urlsplit(reference).path)).resolve()
    if not path.is_relative_to(root.resolve()) or not path.is_file():
        raise ValueError(f"missing or out-of-tree documentation asset: {reference}")
    return path


def data_uri(path):
    kind = mimetypes.guess_type(path.name)[0] or "application/octet-stream"
    return f"data:{kind};base64,"+base64.b64encode(path.read_bytes()).decode("ascii")


def inline_css(path, root):
    css = path.read_text(encoding="utf-8")
    def replace(match):
        reference = match.group(2)
        if reference.startswith(("data:", "#")):
            return match.group(0)
        if urlsplit(reference).scheme or reference.startswith("//"):
            raise ValueError(f"portable CSS cannot depend on remote asset: {reference}")
        # CSS URLs are relative to the stylesheet, not the HTML document.
        target = (path.parent/unquote(urlsplit(reference).path)).resolve()
        if not target.is_relative_to(root.resolve()) or not target.is_file():
            raise ValueError(f"missing CSS asset: {reference}")
        return 'url("'+data_uri(target)+'")'
    return re.sub(r"url\(\s*(['\"]?)(.*?)\1\s*\)", replace, css)


def render_math(soup):
    nodes = soup.select(".math")
    formulas = []
    for node in nodes:
        tex = node.get_text().strip()
        display = tex.startswith(r"\[")
        if tex.startswith((r"\[", r"\(")):
            tex = tex[2:-2]
        formulas.append({"tex": tex, "display": display})
    node_binary = os.environ.get("CODEX_PRIMARY_RUNTIME_NODE") or shutil.which("node")
    if not node_binary:
        raise RuntimeError("Node.js is required to pre-render the portable manual's formulas")
    katex = Path(sphinxcontrib.katex.__file__).parent/"katex.min.js"
    program = """
const fs = require('fs');
const katex = require(process.argv[1]);
try {
  const formulas = JSON.parse(fs.readFileSync(0, 'utf8'));
  const rendered = formulas.map(f => katex.renderToString(f.tex, {
    displayMode: f.display, output: 'mathml', throwOnError: true
  }));
  process.stdout.write(JSON.stringify(rendered));
} catch (error) {
  process.stderr.write(error.message + '\\n');
  process.exitCode = 1;
}
"""
    result = subprocess.run([node_binary, "-e", program, str(katex)],
                            input=json.dumps(formulas), text=True, capture_output=True, check=True)
    for node, markup in zip(nodes, json.loads(result.stdout)):
        replacement = BeautifulSoup(markup, "html.parser")
        node.clear()
        node.append(replacement)
    return len(nodes)


def validate(soup):
    """Check actual delivered navigation/resources, not just build success."""
    identifiers = {node["id"] for node in soup.select("[id]")}
    anchors = 0
    for link in soup.find_all("a", href=True):
        reference = link["href"]
        if reference.startswith("#"):
            if unquote(reference[1:]) not in identifiers:
                raise ValueError(f"broken manual anchor: {reference}")
            anchors += 1
        elif not urlsplit(reference).scheme and not reference.startswith("//"):
            raise ValueError(f"manual depends on another local page: {reference}")
    if soup.find("script") or soup.find("form"):
        raise ValueError("portable documentation must not require scripts or an external search form")
    for node in soup.find_all(src=True):
        if not node["src"].startswith("data:"):
            raise ValueError(f"manual depends on another resource: {node['src']}")
    for node in soup.find_all("link"):
        if node.get("rel") != ["canonical"]:
            raise ValueError("portable manual contains an external link resource")
    languages = soup.select(".language-switch a")
    if [a.get_text(strip=True) for a in languages] != ["English", "Čeština"]:
        raise ValueError("missing language navigation")
    if languages[1].get("lang") != "cs" or languages[1].get("href") != "#document-cesky":
        raise ValueError("Czech language link must target the embedded Czech guide")
    if not soup.find("math"):
        raise ValueError("manual lacks pre-rendered native MathML")
    return anchors


def make_portable(source, destination):
    source, destination = Path(source), Path(destination)
    soup = BeautifulSoup(source.read_text(encoding="utf-8"), "html.parser")
    root = source.parent.resolve()
    soup.body["id"] = "manual-top"
    equations = render_math(soup)
    for script in soup.find_all("script"):
        script.decompose()
    for search in soup.select("search, form, #searchbox"):
        search.decompose()
    for link in list(soup.find_all("link")):
        if "stylesheet" in link.get("rel", []):
            style = soup.new_tag("style")
            style.string = inline_css(asset(root, link["href"]), root)
            link.replace_with(style)
        else:
            link.decompose()
    for image in soup.find_all("img", src=True):
        image["src"] = data_uri(asset(root, image["src"]))
    # A singlehtml template's pathto() collapses both language paths to '#'.
    # Give these controls explicit chapter targets in the assembled manual.
    for link in soup.select('.language-switch a[lang="en"]'):
        link["href"] = "#manual-top"
        link.attrs.pop("aria-current", None)
    for link in soup.select('.language-switch a[lang="cs"]'):
        link["href"] = "#document-cesky"
    czech = soup.find(id="document-cesky")
    if czech is not None:
        czech["lang"] = "cs"
    identifiers = {node["id"] for node in soup.select("[id]")}
    for link in list(soup.find_all("a", href=True)):
        url = urlsplit(link["href"])
        if url.scheme or link["href"].startswith("//"):
            continue
        if url.path in ("", source.name):
            fragment = url.fragment
            if fragment in ("", "document-index"):
                fragment = "manual-top"
            # Sphinx cross-document API references can retain two fragments
            # in a singlehtml build: #document-api#qualified.member. The
            # actual merged document gives the member its own global id.
            if fragment not in identifiers and "#" in fragment:
                member = fragment.rsplit("#", 1)[-1]
                if member in identifiers:
                    fragment = member
            link["href"] = "#"+fragment
        elif url.path.startswith(("_modules/", "_sources/")):
            link.decompose()
        elif url.path in ("genindex.html", "search.html", "py-modindex.html"):
            link.decompose()
        else:
            raise ValueError(f"unexpected separate page in single-document build: {url.path}")
    # A source-code [source] link may already be represented as an inner span.
    for link in soup.select(".viewcode-link"):
        link.decompose()
    anchors = validate(soup)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(str(soup), encoding="utf-8")
    print(f"Portable manual: {destination} ({anchors} checked anchors, {equations} pre-rendered equations)")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT/"docs/_build/padic-tools-manual.html")
    arguments = parser.parse_args()
    build = ROOT/"docs/_build/singlehtml"
    subprocess.run([sys.executable, "-m", "sphinx", "-W", "--keep-going", "-b", "singlehtml",
                    str(ROOT/"docs"), str(build)], check=True)
    make_portable(build/"index.html", arguments.output)


if __name__ == "__main__":
    main()
