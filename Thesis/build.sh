#!/bin/bash

set -e  # stop on error

echo "Generating meta.tex from config.json..."
python3 generate_meta.py config.json

run_process () {
    echo "> $1"
    eval "$1"
}

if command -v latexmk >/dev/null 2>&1; then
    echo "Using latexmk to build PDF..."
    run_process "latexmk -pdf main.tex"
    echo "Build finished: main.pdf"

    echo "running live run"
    run_process "latexmk -pdf -pvc main.tex"
else
    echo "latexmk not found, falling back to pdflatex + biber..."
    run_process "pdflatex -interaction=nonstopmode main.tex"
    run_process "biber main"
    run_process "pdflatex -interaction=nonstopmode main.tex"
    run_process "pdflatex -interaction=nonstopmode main.tex"
fi