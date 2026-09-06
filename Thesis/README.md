# LaTeX Paper Template

This project contains a configurable LaTeX paper scaffold.

Quick start (Windows PowerShell):

1. Edit `config.json` to set title, author, document class, and packages.
2. Run the metadata generator (creates `meta.tex`):

```powershell
python .\generate_meta.py config.json
```

3. Build the PDF (uses `latexmk` if available, otherwise pdflatex + biber/bibtex):

```powershell
.\build.ps1
```

- For Linux Os: 
```powershell
chmod +x build.sh
# then: 
./build.sh
```
- Make sure you have latexmk and biber and other files installes:

```powershell
sudo apt update
# then
sudo apt install -y latexmk texlive-latex-extra biber
```
- remember that now everything is case- sensitive , make sure to check for exact formats, capital cases, and the like.
Files of interest:
- `config.json`: project metadata and packages
- `generate_meta.py`: turns `config.json` into `meta.tex`
- `main.tex`: main document (includes `meta.tex`)
- `sections/`: example sections
- `bibliography.bib`: sample bib file
- `build.ps1`: Windows build script



### Latex language cheetsheet
# LaTeX Comprehensive Guide

## 1. Document Structure

### Basic Document
```latex
\documentclass[options]{class}
\usepackage{package1, package2}

\begin{document}
% Content here
\end{document}
```

### Common Document Classes
```latex
\documentclass{article}    % Articles
\documentclass{report}     % Multi-chapter reports
\documentclass{book}       % Books
\documentclass{beamer}     % Presentations
\documentclass{letter}     % Letters
```

### Options
```latex
\documentclass[11pt, a4paper, twocolumn]{article}
% 10pt, 11pt, 12pt - font sizes
% letterpaper, a4paper - paper size
% twocolumn - two columns
% landscape - landscape orientation
% titlepage - separate title page
```

## 2. Text Formatting

### Font Styles
```latex
\textbf{Bold Text}
\textit{Italic Text}
\textsl{Slanted Text}
\textsc{Small Caps}
\underline{Underlined}
\emph{Emphasized} % Italic by default, bold in bold context
```

### Font Sizes
```latex
{\tiny tiny text}
{\scriptsize scriptsize text}
{\footnotesize footnotesize text}
{\small small text}
{\normalsize normal text}
{\large large text}
{\Large Larger text}
{\LARGE LARGE text}
{\huge huge text}
{\Huge Huge text}
```

### Font Families
```latex
\textrm{Roman} \textsf{Sans Serif} \texttt{Typewriter}
\rmfamily Roman \sffamily Sans Serif \ttfamily Typewriter
```

### Colors (requires `xcolor` package)
```latex
\usepackage{xcolor}
\textcolor{red}{Red text}
\colorbox{yellow}{Highlighted}
\fcolorbox{blue}{gray}{Bordered box}
```

## 3. Mathematics

### Inline Math
```latex
$E = mc^2$
```

### Display Math
```latex
\[ E = mc^2 \]
\begin{equation}
E = mc^2
\end{equation}
```

### Aligned Equations
```latex
\begin{align}
    y &= mx + b \\
    f(x) &= x^2 + 2x + 1
\end{align}
```

### Mathematical Symbols

#### Greek Letters
```latex
\alpha \beta \gamma \Gamma \delta \Delta \pi \Pi \omega \Omega
```

#### Operators
```latex
\times \div \pm \mp \cdot \ast \bullet \sum \prod \int \oint
\lim \max \min \sup \inf \det \log \ln \exp \sin \cos \tan
```

#### Relations
```latex
= \neq < > \leq \leq \geq \geq \approx \sim \propto \equiv
\in \subset \subseteq \supset \supseteq \cup \cap
```

#### Arrows
```latex
\to \rightarrow \Rightarrow \leftrightarrow \Leftrightarrow
\uparrow \downarrow \mapsto \longmapsto
```

#### Special Symbols
```latex
\infty \forall \exists \emptyset \nabla \partial \hbar
\ell \Re \Im \aleph \dagger \ddagger \S \P
```

#### Fractions and Roots
```latex
\frac{a}{b} \tfrac{a}{b} \dfrac{a}{b}
\sqrt{x} \sqrt[n]{x}
```

#### Matrices
```latex
\begin{matrix} a & b \\ c & d \end{matrix}
\begin{pmatrix} a & b \\ c & d \end{pmatrix}
\begin{bmatrix} a & b \\ c & d \end{bmatrix}
\begin{vmatrix} a & b \\ c & d \end{vmatrix}
```

## 4. Lists

### Bulleted List
```latex
\begin{itemize}
    \item First item
    \item Second item
    \begin{itemize}
        \item Nested item
    \end{itemize}
\end{itemize}
```

### Numbered List
```latex
\begin{enumerate}
    \item First item
    \item Second item
    \begin{enumerate}
        \item Nested item
    \end{enumerate}
\end{enumerate}
```

### Description List
```latex
\begin{description}
    \item[Term 1] Description 1
    \item[Term 2] Description 2
\end{description}
```

## 5. Tables

### Basic Table
```latex
\begin{tabular}{|l|c|r|}
    \hline
    Left & Center & Right \\
    \hline
    A & B & C \\
    \hline
\end{tabular}
```

### Table with Caption
```latex
\begin{table}[h]
    \centering
    \caption{Table Title}
    \begin{tabular}{lll}
        \toprule
        Header 1 & Header 2 & Header 3 \\
        \midrule
        Data 1 & Data 2 & Data 3 \\
        Data 4 & Data 5 & Data 6 \\
        \bottomrule
    \end{tabular}
    \label{tab:my_table}
\end{table}
```

### Column Alignment
- `l` - left aligned
- `c` - centered
- `r` - right aligned
- `p{3cm}` - fixed width paragraph
- `|` - vertical line
- `\hline` - horizontal line

## 6. Figures

### Including Images
```latex
\usepackage{graphicx}
\begin{figure}[h]
    \centering
    \includegraphics[width=0.5\textwidth]{image.jpg}
    \caption{Image caption}
    \label{fig:my_image}
\end{figure}
```

### Image Options
```latex
\includegraphics[width=0.8\textwidth]{file}
\includegraphics[height=3cm]{file}
\includegraphics[scale=0.5]{file}
\includegraphics[angle=45]{file}
\includegraphics[trim=left bottom right top]{file}
```

## 7. References and Citations

### Labels and References
```latex
\section{Introduction}\label{sec:intro}
See Section~\ref{sec:intro} for details.
Equation~\eqref{eq:emc} shows the relationship.
```

### Bibliography with BibTeX
```latex
\cite{key}          % Basic citation
\citep{key}         % Parenthetical citation
\citet{key}         % Textual citation
\citeauthor{key}    % Author only
\citeyear{key}      % Year only
```

## 8. Document Elements

### Sections
```latex
\section{Section Title}
\subsection{Subsection Title}
\subsubsection{Subsubsection Title}
\paragraph{Paragraph Title}
\subparagraph{Subparagraph Title}
```

### Footnotes
```latex
Main text\footnote{Footnote text here.}
```

### Verbatim Text
```latex
\begin{verbatim}
Code or text exactly as typed
\end{verbatim}
```

### Quotes and Quotations
```latex
\begin{quote}
    Short quotation.
\end{quote}

\begin{quotation}
    Longer quotation with indentation.
\end{quotation}
```

## 9. Special Characters and Spacing

### Escaping Special Characters
```latex
\% \$ \& \# \_ \{ \} \~{} \^{} \textbackslash
```

### Dashes
```latex
-   % hyphen (short dash)
--  % en-dash (for ranges)
--- % em-dash (for punctuation)
```

### Spacing
```latex
\,    % thin space
\:    % medium space
\;    % thick space
\quad % 1 em space
\qquad % 2 em space
\!    % negative thin space
```

### Line Breaks
```latex
\\    % line break
\\[10pt] % line break with space
\newline % new line
\pagebreak % page break
\clearpage % clear page
```

## 10. Custom Commands and Environments

### New Commands
```latex
\newcommand{\R}{\mathbb{R}}
\newcommand{\abs}[1]{\left|#1\right|}
\newcommand{\pderiv}[2]{\frac{\partial #1}{\partial #2}}
```

### Renew Commands
```latex
\renewcommand{\vec}[1]{\mathbf{#1}}
```

### New Environments
```latex
\newenvironment{myenv}
    {\begin{center}\itshape}
    {\end{center}}
```

## 11. Common Packages

### Essential Packages
```latex
\usepackage{amsmath}        % Advanced math
\usepackage{amssymb}        % Additional symbols
\usepackage{graphicx}       % Graphics
\usepackage{hyperref}       % Hyperlinks
\usepackage{geometry}       % Page layout
\usepackage{fancyhdr}       % Headers and footers
\usepackage{booktabs}       % Professional tables
\usepackage{array}          % Enhanced tables
\usepackage{xcolor}         % Colors
\usepackage{listings}       % Code listings
\usepackage{biblatex}       % Bibliography
\usepackage{float}          % Better float control
\usepackage{siunitx}        % SI units
```

### Hyperref Setup
```latex
\usepackage{hyperref}
\hypersetup{
    colorlinks=true,
    linkcolor=blue,
    filecolor=magenta,      
    urlcolor=cyan,
    citecolor=green,
    pdftitle={Your Title},
    pdfauthor={Your Name}
}
```

## 12. Troubleshooting Common Issues

### "Missing $ inserted"
- Use math mode for math symbols: `$E = mc^2$`
- Escape underscores in text: `variable\_name`

### "Undefined control sequence"
- Check package is loaded
- Check spelling of command
- Add missing package

### "File not found"
- Check file paths are correct
- Use relative paths
- Ensure file extensions are included

### "Overfull \hbox"
- Line too long
- Use hyphenation or manual line breaks
- Adjust margins with `geometry` package

### Float Placement Issues
```latex
\begin{figure}[h!]   % h! forces here if possible
\begin{figure}[H]    % Requires float package, forces here
```

## 13. Tips and Best Practices

1. **Use consistent naming** for labels: `fig:`, `tab:`, `eq:`, `sec:`
2. **Keep .tex files organized** in folders
3. **Comment your code** with `%`
4. **Compile frequently** to catch errors early
5. **Use version control** (Git) for LaTeX projects
6. **Separate content** into multiple files using `\input` or `\include`
7. **Use BibTeX** for references management
8. **Check spelling** and grammar
9. **Use `\usepackage[utf8]{inputenc}`** for UTF-8 support
10. **Include a Makefile** or compile script for complex projects

## 14. Quick Reference Card

```latex
% Comments
% Math mode: $...$ or \[...\]
% Bold: \textbf{...}
% Italic: \textit{...}
% Sections: \section{...}
% References: \label{...} and \ref{...}
% Citations: \cite{...}
% Figures: \includegraphics[...]{...}
% Tables: \begin{tabular}{...}
% Lists: \begin{itemize}, \begin{enumerate}
% Footnotes: \footnote{...}
% URLs: \url{...} or \href{...}{...}
```

This guide covers the most essential LaTeX commands and formatting options. For more detailed information, consult the comprehensive documentation for specific packages or the LaTeX Wikibook.