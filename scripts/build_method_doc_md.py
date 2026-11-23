"""Convert the generated Word report to GitHub-flavoured Markdown.

The .docx is the source of truth: this walks the body in document order,
keeps headings/captions/bullets/quotes, and turns every real table into a
GFM pipe table.  The abstract box (a 1x1 table) becomes a block quote.

Usage::

    python scripts/build_method_doc_md.py                 # rebuild docs/...md
    python scripts/build_method_doc_md.py in.docx out.md  # explicit paths
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
from docx.table import Table
from docx.text.paragraph import Paragraph

ROOT = Path(__file__).resolve().parents[1]
SRC = (
    Path(sys.argv[1])
    if len(sys.argv) > 1
    else ROOT / "docs" / "AutoSMJ_DFT方法说明_v1.2.docx"
)
DST = (
    Path(sys.argv[2])
    if len(sys.argv) > 2
    else ROOT / "docs" / "AutoSMJ_DFT方法说明_v1.2.md"
)


def esc(text: str) -> str:
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def rich(paragraph: Paragraph) -> str:
    """Paragraph text with **bold** segments preserved."""
    segments: list[list] = []
    for run in paragraph.runs:
        text = esc(run.text)
        if not text:
            continue
        bold = bool(run.bold)
        if segments and segments[-1][1] == bold:
            segments[-1][0] += text
        else:
            segments.append([text, bold])

    out = []
    for text, bold in segments:
        if bold and text.strip():
            out.append(f"**{text}**")
        else:
            out.append(text)
    s = "".join(out)
    # A bold lead ending in a full-width colon should breathe before the text.
    s = re.sub(r"\*\*([^*]*：)\*\*(?=\S)", r"**\1** ", s)
    return s


def cell_text(cell) -> str:
    parts = []
    for p in cell.paragraphs:
        t = p.text.strip()
        if t:
            parts.append(esc(t))
    return "<br>".join(parts).replace("|", r"\|")


def column_alignment(table: Table, j: int) -> str:
    for row in table.rows[1:]:
        if j >= len(row.cells):
            continue
        p = row.cells[j].paragraphs[0]
        if p.alignment == WD_ALIGN_PARAGRAPH.RIGHT:
            return "---:"
        if p.alignment == WD_ALIGN_PARAGRAPH.CENTER:
            return ":---:"
        if p.alignment == WD_ALIGN_PARAGRAPH.LEFT:
            return "---"
    return "---"


def md_table(table: Table) -> str:
    rows = []
    for row in table.rows:
        rows.append([cell_text(c) for c in row.cells])
    n = max(len(r) for r in rows)
    rows = [r + [""] * (n - len(r)) for r in rows]

    header = rows[0]
    aligns = [column_alignment(table, j) for j in range(n)]
    lines = [
        "| " + " | ".join(header) + " |",
        "| " + " | ".join(aligns) + " |",
    ]
    for row in rows[1:]:
        lines.append("| " + " | ".join(row) + " |")
    return "\n".join(lines)


def abstract_block(table: Table) -> str:
    blocks: list[str] = []
    for p in table.rows[0].cells[0].paragraphs:
        t = p.text.strip()
        if not t:
            continue
        if t == "摘要":
            blocks.append("**摘要**")
        else:
            blocks.append(rich(p))
    return "\n>\n".join("> " + b if b else ">" for b in blocks)


doc = Document(SRC)
items: list[tuple[str, object]] = []
for child in doc.element.body.iterchildren():
    if child.tag == qn("w:p"):
        items.append(("p", Paragraph(child, doc)))
    elif child.tag == qn("w:tbl"):
        items.append(("t", Table(child, doc)))

# ---- cover ----------------------------------------------------------------
first_h = next(
    i for i, (kind, obj) in enumerate(items)
    if kind == "p" and obj.style.name == "Heading 1"
)
cover_ps = [o for k, o in items[:first_h] if k == "p" and o.text.strip()]
abstract_tables = [o for k, o in items[:first_h] if k == "t"]

title = cover_ps[1].text.replace("\n", " ").strip()
kicker = cover_ps[0].text.strip()
subtitle = cover_ps[2].text.strip()
meta_lines = []
for p in cover_ps[3:7]:
    label, _, value = p.text.partition("\u3000")
    meta_lines.append(f"- **{label.strip()}**：{value.strip()}")
disclaimer = cover_ps[7].text.strip()

out: list[str] = []
out.append(f"# {title}")
out.append(f"*{kicker}*")
out.append(f"**{subtitle}**")
out.append("\n".join(meta_lines))
out.append(f"*{disclaimer}*")
for table in abstract_tables:
    out.append(abstract_block(table))

# ---- body ------------------------------------------------------------------
for kind, obj in items[first_h:]:
    if kind == "t":
        table: Table = obj  # type: ignore[assignment]
        if len(table.rows) == 1 and len(table.columns) == 1 and \
                table.rows[0].cells[0].text.strip().startswith("摘要"):
            out.append(abstract_block(table))
        else:
            out.append(md_table(table))
        continue

    p: Paragraph = obj  # type: ignore[assignment]
    text = p.text.strip()
    if not text:
        continue
    if p.style.name == "Heading 1":
        out.append("## " + re.sub(r"\s+", " ", text))
        continue
    if p.style.name == "Heading 2":
        out.append("### " + re.sub(r"\s+", " ", text))
        continue
    if text.startswith("表 ") and p.alignment == WD_ALIGN_PARAGRAPH.CENTER:
        out.append("**" + re.sub(r"\s+", " ", text) + "**")
        continue
    if text.startswith("▪"):
        body = rich(p).lstrip("▪").lstrip()
        out.append("- " + body)
        continue
    ppr = p._p.pPr
    shd = ppr.find(qn("w:shd")) if ppr is not None else None
    fill = shd.get(qn("w:fill")) if shd is not None else None
    if fill == "F2F4F7":
        out.append("> " + rich(p))
        continue
    out.append(rich(p))

rendered = ""
for i, block in enumerate(out):
    if i == 0:
        rendered = block
        continue
    sep = "\n" if (block.startswith("- ") and out[i - 1].startswith("- ")) else "\n\n"
    rendered += sep + block

DST.write_bytes((rendered + "\n").encode("utf-8"))
print("saved:", DST)
print("chars:", len(DST.read_text(encoding="utf-8")))
