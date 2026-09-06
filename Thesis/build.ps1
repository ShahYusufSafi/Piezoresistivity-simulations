Write-Host "Generating meta.tex from config.json..."
python .\generate_meta.py config.json

function Run-Process($cmd) {
    Write-Host "> $cmd"
    & cmd /c $cmd
    if ($LASTEXITCODE -ne 0) { throw "Command failed: $cmd" }
}

if (Get-Command latexmk -ErrorAction SilentlyContinue) {
    Write-Host "Using latexmk to build PDF..."
    Run-Process 'latexmk -pdf main.tex'
    Write-Host "Build finished: main.pdf"
}
else {
    Write-Host "latexmk not found, falling back to pdflatex + biber..."
    Run-Process "pdflatex -interaction=nonstopmode main.tex"
    Run-Process "biber main"
    Run-Process "pdflatex -interaction=nonstopmode main.tex"
    Run-Process "pdflatex -interaction=nonstopmode main.tex"
}
