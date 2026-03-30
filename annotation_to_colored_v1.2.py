#!/usr/bin/env python3
#This script has been created with the help of ChatGPT by Lennart Johansson 28-03-2026
import csv
import argparse

# kleur mapping
COLOR_MAP = {
    "A": "#4CAF50",  # green
    "C": "#2196F3",  # blue
    "G": "#FFEB3B",  # yellow
    "T": "#F44336",  # red
}

INTERRUPTION_COLOR = "#BDBDBD"  # grey


def color_sequence(seq, motif):
    result = []
    i = 0
    mlen = len(motif)

    while i < len(seq):
        base = seq[i]

        if base.isupper():
            # strict match
            if seq[i:i+mlen] == motif:
                for j in range(mlen):
                    b = seq[i+j]
                    color = COLOR_MAP.get(b, "#FFFFFF")
                    result.append(
                        f'<span style="background-color:{color}; padding:2px; font-family:monospace;">{b}</span>'
                    )
                i += mlen
            else:
                # interruption → grey per base
                result.append(
                    f'<span style="background-color:{INTERRUPTION_COLOR}; padding:2px; font-family:monospace;">{base}</span>'
                )
                i += 1
        else:
            # flank (lowercase)
            result.append(base)
            i += 1

    return "".join(result)

def safe_repeat_bp(x):
    try:
        return int(x["repeat_bp"])
    except:
        return -1  # ensures NA rows go to bottom


def main():
    parser = argparse.ArgumentParser(
        description="Convert annotation.tsv to colored HTML with strict repeat highlighting"
    )
    parser.add_argument("input", help="Input TSV file (annotation.tsv)")
    parser.add_argument("output", help="Output HTML file (annotation.colored.html)")
    args = parser.parse_args()

    # inlezen
    with open(args.input, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
        fieldnames_original = reader.fieldnames

    # sorteren op repeat_bp (langste boven)
    rows.sort(key=safe_repeat_bp, reverse=True)

    # seq_between_coords naar voren halen
    fieldnames = ["seq_between_coords"] + [
        f for f in fieldnames_original if f != "seq_between_coords"
    ]

    # schrijven
    with open(args.output, "w") as f:
        f.write("""
<html>
<head>
<style>
table {
    border-collapse: collapse;
    table-layout: auto;
    font-family: monospace;
}

td, th {
    border: 1px solid black;
    padding: 4px;
    white-space: nowrap;
}

th {
    position: sticky;
    top: 0;
    background: white;
}
</style>
</head>
<body>
<div style="overflow-x:auto;">
<table>
""")

        # header
        f.write("<tr>")
        for col in fieldnames:
            f.write(f"<th>{col}</th>")
        f.write("</tr>\n")

        # rows
        for row in rows:
            f.write("<tr>")
            motif = row["target_repeat"]

            for col in fieldnames:
                val = row[col]

                if col == "seq_between_coords":
                    if val and val != "NA":
                        val = color_sequence(val, motif)
                    else:
                        val = ""  # empty block

                elif col == "read_id":
                    val = f'<span style="font-size:10px;">{val}</span>'

                f.write(f"<td>{val}</td>")

            f.write("</tr>\n")

        f.write("""
</table>
</div>
</body>
</html>
""")

if __name__ == "__main__":
    main()
