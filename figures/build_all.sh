#!/usr/bin/env bash
# Build every figure standalone.  ./build_all.sh [FigN ...]
#
# Each FigN/main.tex reproduces the text block the paper gives that figure, so
# the PDF it produces is the same size as in the article. See README.md.
set -u
cd "$(dirname "$0")"
targets=("$@"); [ ${#targets[@]} -eq 0 ] && targets=(Fig1 Fig2 Fig3 Fig4 Fig5 Fig6 Fig7 Fig8 Fig9 Fig10)
fail=0
for d in "${targets[@]}"; do
    [ -f "$d/main.tex" ] || { printf '  %-6s no main.tex\n' "$d"; continue; }
    if ( cd "$d" && pdflatex -interaction=nonstopmode -halt-on-error main.tex >build.log 2>&1 ); then
        u=$(grep -c 'Undefined control sequence' "$d/build.log" || true)
        printf '  %-6s ok   main.pdf' "$d"
        [ "$u" -gt 0 ] && printf '   WARNING: %s undefined control sequence(s)' "$u"
        printf '\n'
    else
        printf '  %-6s FAILED (see %s/build.log)\n' "$d" "$d"; fail=1
    fi
done
exit $fail
