"""Assemble the single-file page: src/page.html + src/model.js + src/swan.js + data/cobs_lightcurves.txt + src/ui.js.
Writes index.html (full document, for GitHub Pages) and dist/artifact.html (body-only form)."""
from pathlib import Path
here = Path(__file__).parent
src = here / "src"
page = (src / "page.html").read_text()
model = (src / "model.js").read_text()
# COBS light curves (data/cobs_lightcurves.txt, CC BY-NC-SA 4.0) as a JS object of preset texts
cobs, key = {}, None
for line in (here / "data" / "cobs_lightcurves.txt").read_text().splitlines():
    if line.startswith("#"): key = line[1:].strip(); cobs[key] = []
    elif line.strip() and not line.startswith("%"): cobs[key].append(line.strip())
import json
code = model + "\n" + (src / "swan.js").read_text() + "\nconst COBS = " + json.dumps({k: "\n".join(v) for k, v in cobs.items()}) + ";\n"
page = page.replace("/*MODEL*/", code).replace("/*UI*/", (src / "ui.js").read_text())
(here / "dist").mkdir(exist_ok=True)
(here / "dist" / "artifact.html").write_text(page)
head, body = page.split('<div class="wrap">', 1)
head = head.replace('<meta name="claude:viewport-fit-sides" content="cover">\n', "")
doc = ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
       '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n'
       + head + '<style>body{margin:0}</style>\n</head>\n<body>\n<div class="wrap">' + body + "\n</body>\n</html>\n")
(here / "index.html").write_text(doc)
print("index.html", len(doc), "bytes; dist/artifact.html", len(page), "bytes")
