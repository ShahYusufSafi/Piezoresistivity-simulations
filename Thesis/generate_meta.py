#!/usr/bin/env python3
"""
Simple generator to turn `config.json` into `meta.tex`.

Usage: python generate_meta.py config.json
"""
import json
import sys
from pathlib import Path


def tex_escape(s: str) -> str:
    return s


def main():
    if len(sys.argv) < 2:
        print("Usage: python generate_meta.py config.json")
        sys.exit(1)
    cfg_path = Path(sys.argv[1])
    cfg = json.loads(cfg_path.read_text(encoding="utf8"))

    out = ["% !TEX root = main.tex"]
    #docclass = cfg.get("documentclass", "article")
    #class_opts = cfg.get("class_options", "11pt")
    #out.append(f"\\documentclass[{class_opts}]{{{docclass}}}")

    # geometry
    geom = cfg.get("geometry_options")
    if geom:
        out.append(f"\\usepackage[{geom}]{{geometry}}")

    # packages
    for pkg in cfg.get("packages", []):
        out.append(f"\\usepackage{{{pkg}}}")

    # bibliography (biblatex)
    bib = cfg.get("bibliography")
    if bib:
        out.append("\\usepackage[backend=biber,style=numeric]{biblatex}")
        out.append(f"\\addbibresource{{{bib}}}")

    # title/author/date
    title = cfg.get("title", "")
    author = cfg.get("author", "")
    date = cfg.get("date", "\\today")

    out.append(f"\\title{{{tex_escape(title)}}}")
    out.append(f"\\author{{{tex_escape(author)}}}")
    out.append(f"\\date{{{date}}}")

    meta_path = Path("meta.tex")
    meta_path.write_text("\n".join(out) + "\n", encoding="utf8")
    print(f"Wrote {meta_path.resolve()}")

if __name__ == "__main__":
    main()
