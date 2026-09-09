"""Optional Word-document report showing original vs. optimized sequences
with per-nucleotide differences highlighted.
"""

import logging
from os import path

try:
    from docx import Document
    from docx.enum.text import WD_COLOR_INDEX
    DOCX_AVAILABLE = True
except ImportError:
    DOCX_AVAILABLE = False

logger = logging.getLogger('eso.report')


def _add_highlighted_runs(paragraph, seq, other_seq):
    """Add `seq` to `paragraph` as a run per contiguous stretch of positions
    that all differ from `other_seq`, or all match it - not one run per
    character. A per-character run is what this used to do: for even a
    moderately long real sequence (hundreds to thousands of nt) that's
    thousands of separate `<w:r>` XML elements in the generated .docx - real,
    measurable overhead and file bloat, not just a style choice, for
    something python-docx already lets a single `add_run` call cover.
    A position past the end of `other_seq` is treated as "not differing"
    (unhighlighted) - matches this module's pre-existing behavior; eso's own
    optimization is always substitution-based (DNAChisel never changes
    sequence length), so `seq`/`other_seq` are the same length in every real
    call from eso.pipeline, and this only matters for a hand-built call with
    mismatched lengths.
    """
    if not seq:
        return

    def differs(i):
        return i < len(other_seq) and seq[i] != other_seq[i]

    start = 0
    current_state = differs(0)
    for i in range(1, len(seq) + 1):
        state = differs(i) if i < len(seq) else not current_state
        if state != current_state:
            run = paragraph.add_run(seq[start:i])
            if current_state:
                run.font.highlight_color = WD_COLOR_INDEX.YELLOW
            start = i
            current_state = state


def create_word_document_with_highlighted_differences(sequences_data, output_path):
    """
    Parameters
    ----------
    sequences_data: list of (sequence_name, original_seq, final_seq) tuples.
    output_path: directory in which to save 'sequence_comparison.docx'.
    """
    if not DOCX_AVAILABLE:
        logger.info(
            "python-docx not available (install the 'docx-report' extra), skipping Word document generation")
        return

    doc = Document()
    doc.add_heading('Sequence Optimization Results', 0)

    for index, (seq_name, original_seq, final_seq) in enumerate(sequences_data):
        doc.add_heading(f'Sequence: {seq_name}', level=1)

        doc.add_heading('Original Sequence:', level=2)
        original_paragraph = doc.add_paragraph()

        doc.add_heading('Final Sequence:', level=2)
        final_paragraph = doc.add_paragraph()

        _add_highlighted_runs(original_paragraph, original_seq, final_seq)
        _add_highlighted_runs(final_paragraph, final_seq, original_seq)

        # was `if seq_name != sequences_data[-1][0]` - compared by name, so a
        # sequence whose name happened to match the true last entry's name
        # (e.g. two files/records that coincidentally share a stem) would
        # wrongly skip its own page break. Compare by position instead.
        if index != len(sequences_data) - 1:
            doc.add_page_break()

    doc_path = path.join(output_path, 'sequence_comparison.docx')
    doc.save(doc_path)
    logger.info(f"Word document saved to: {doc_path}")
