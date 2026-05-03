"""``adataviz.pl`` -- entry point for all plot functions.

Real implementations live in topic-specific modules:

==============================  =========================================
Module                          Plots
==============================  =========================================
``scatter``                     plot_categorical / plot_gene
``expression``                  boxplot / stacked_violinplot / gene_dotplot
``pseudotime``                  plot_pseudotime
``composition``                 rose / ring / pie / area / dot / trend / bar
``flow``                        sankey / chord
``sets``                        venn / upset
``complex_heatmap``             complex_heatmap / complex_dotplot
``interactive_stats``           plotly versions of all of the above
==============================  =========================================
"""

from __future__ import annotations

from ._utils import use_scientific_style, show_fig
from .scatter import plot_categorical, plot_gene, assign_default_colors
from .expression import (
    boxplot,
    stacked_violinplot,
    gene_dotplot,
    get_genes_mean_frac,
)
from .pseudotime import plot_pseudotime
from .composition import (
    rose_plot,
    ring_plot,
    pie_plot,
    area_plot,
    dot_plot,
    trend_plot,
    bar_plot,
)
from .flow import sankey_plot, chord_plot
from .sets import venn_plot, upset_plot
from .complex_heatmap import complex_heatmap, complex_dotplot
from .interactive_stats import (
    interactive_rose,
    interactive_ring,
    interactive_pie,
    interactive_area,
    interactive_trend,
    interactive_sankey,
    interactive_embedding,
    interactive_categorical,
    interactive_gene,
    interactive_dotHeatmap,
    interactive_boxplot,
    interactive_violin,
    interactive_stacked_bar,
    interactive_dot_plot,
    interactive_chord,
    interactive_upset,
    interactive_complex_heatmap,
    interactive_complex_dotplot,
)

__all__ = [
    "use_scientific_style",
    "show_fig",
    "plot_categorical",
    "plot_gene",
    "assign_default_colors",
    "boxplot",
    "stacked_violinplot",
    "gene_dotplot",
    "get_genes_mean_frac",
    "plot_pseudotime",
    "rose_plot",
    "ring_plot",
    "pie_plot",
    "area_plot",
    "dot_plot",
    "trend_plot",
    "bar_plot",
    "sankey_plot",
    "chord_plot",
    "venn_plot",
    "upset_plot",
    "complex_heatmap",
    "complex_dotplot",
    "interactive_rose",
    "interactive_ring",
    "interactive_pie",
    "interactive_area",
    "interactive_trend",
    "interactive_sankey",
    "interactive_embedding",
    "interactive_categorical",
    "interactive_gene",
    "interactive_dotHeatmap",
    "interactive_boxplot",
    "interactive_violin",
    "interactive_stacked_bar",
    "interactive_dot_plot",
    "interactive_chord",
    "interactive_upset",
    "interactive_complex_heatmap",
    "interactive_complex_dotplot",
]
