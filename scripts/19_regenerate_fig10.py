"""Regenerate figures/fig10_caterpillar_maximizers.{png,pdf} with explicit
ASCII hyphen-minus in the title so the W(T; -1, -1, 1) parameter triple is
unambiguously visible in the compiled MATCH-template PDF at A5 / 10pt.

The earlier hand-made figure used Unicode minus U+2212, which compiles fine
in the PNG but renders extremely thin in the embedded MATCH PDF at A5 font
size, making the title look like W(T; 1, 1, 1).

Trees plotted (the three rooted maximizers for n = 10, 11, 12):
  n = 10: balanced double-spider DS(3, 3),    graph6 = Ip_I?D??G
  n = 11: spider S^(2)_5  (hub degree 5),     graph6 = JqD?I?@O??_
  n = 12: balanced double-spider DS(4, 3),    graph6 = Kp_I?D??I??@

We build the trees by direct construction (hub + length-2 spokes / double-hub
+ length-2 spokes) rather than from graph6 strings, so the layout is
deterministically hierarchical and the rendering is stable across matplotlib
versions.
"""
from __future__ import annotations
import os, sys
import matplotlib.pyplot as plt
import matplotlib
import networkx as nx

# Force ASCII hyphen-minus and a serif math font that matches the body text.
matplotlib.rcParams.update({
    "font.family": "serif",
    "axes.unicode_minus": False,     # critical: use ASCII hyphen-minus, not U+2212
    "mathtext.fontset": "cm",
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
})

PROJECT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
FIG_DIR = os.path.join(PROJECT, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# ----- Tree constructors (rooted layouts) -----

def spider(hub_degree: int) -> nx.Graph:
    """Spider S^(2)_h: one hub of degree h with h length-2 spokes."""
    G = nx.Graph()
    G.add_node(0)              # hub
    next_id = 1
    for _ in range(hub_degree):
        mid = next_id
        leaf = next_id + 1
        G.add_edges_from([(0, mid), (mid, leaf)])
        next_id += 2
    return G


def double_spider(a: int, b: int) -> nx.Graph:
    """Double-spider DS(a, b): two hubs of degrees a, b joined by an edge,
    each carrying length-2 spokes to leaves."""
    G = nx.Graph()
    hub_a, hub_b = 0, 1
    G.add_edge(hub_a, hub_b)
    next_id = 2
    # Hub a has a-1 length-2 spokes (one bond used by the hub-hub edge)
    for _ in range(a - 1):
        mid = next_id
        leaf = next_id + 1
        G.add_edges_from([(hub_a, mid), (mid, leaf)])
        next_id += 2
    # Hub b has b-1 length-2 spokes
    for _ in range(b - 1):
        mid = next_id
        leaf = next_id + 1
        G.add_edges_from([(hub_b, mid), (mid, leaf)])
        next_id += 2
    return G


# ----- Layouts: tailored to each tree shape, no crossing edges -----

def spider_layout(hub_degree: int) -> dict:
    """Pin the hub at the top; fan out h length-2 spokes symmetrically below.
    Returns positions on the nodes produced by spider()."""
    pos = {0: (0.0, 2.0)}  # hub
    h = hub_degree
    if h == 1:
        xs = [0.0]
    else:
        xs = [-2.5 + i * 5.0 / (h - 1) for i in range(h)]
    next_id = 1
    for x in xs:
        mid = next_id
        leaf = next_id + 1
        pos[mid] = (x, 1.0)
        pos[leaf] = (x, 0.0)
        next_id += 2
    return pos


def double_spider_layout(a: int, b: int) -> dict:
    """Two hubs side-by-side at the top; their length-2 spokes fan out
    *outward* (hub_a left, hub_b right), so no edge ever crosses another.
    Returns positions on the nodes produced by double_spider()."""
    sep = 1.5
    pos = {0: (-sep / 2, 2.0), 1: (sep / 2, 2.0)}  # hub_a and hub_b

    # Hub a's a-1 length-2 spokes fan out to the LEFT
    next_id = 2
    if a >= 2:
        # x positions decrease from -sep/2 leftward
        a_spread = max(2.5, 0.9 * (a - 1))
        # leftmost at -sep/2 - a_spread, rightmost at -sep/2 - 0.5
        if a == 2:
            xs_a = [-sep / 2 - 1.2]
        else:
            xs_a = [-sep / 2 - 0.5 - i * a_spread / (a - 2)
                    for i in range(a - 1)] if a > 2 else [-sep / 2 - 1.2]
        for x in xs_a:
            mid = next_id
            leaf = next_id + 1
            pos[mid] = (x, 1.0)
            pos[leaf] = (x, 0.0)
            next_id += 2

    # Hub b's b-1 length-2 spokes fan out to the RIGHT
    if b >= 2:
        b_spread = max(2.5, 0.9 * (b - 1))
        if b == 2:
            xs_b = [sep / 2 + 1.2]
        else:
            xs_b = [sep / 2 + 0.5 + i * b_spread / (b - 2)
                    for i in range(b - 1)] if b > 2 else [sep / 2 + 1.2]
        for x in xs_b:
            mid = next_id
            leaf = next_id + 1
            pos[mid] = (x, 1.0)
            pos[leaf] = (x, 0.0)
            next_id += 2
    return pos


def draw_tree(ax, G: nx.Graph, pos: dict, title: str):
    # Color: internal vertices green, leaves cream
    leaf_color = "#fff5d6"     # pale cream
    internal_color = "#a8d8b0" # soft green
    node_colors = []
    edge_colors = []
    for n in G.nodes():
        if G.degree(n) == 1:
            node_colors.append(leaf_color)
        else:
            node_colors.append(internal_color)

    nx.draw_networkx_edges(G, pos, ax=ax, edge_color="black", width=1.2)
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=node_colors,
                           node_size=420, edgecolors="black", linewidths=1.0)
    ax.set_title(title, fontsize=10)
    ax.axis("off")
    ax.set_aspect("equal")


# ----- Main figure -----

def main():
    fig, axes = plt.subplots(1, 3, figsize=(13, 5.5))

    # n = 10: DS(3, 3)
    G10 = double_spider(3, 3)
    pos10 = double_spider_layout(3, 3)
    draw_tree(
        axes[0], G10, pos10,
        title=(r"$n = 10$, degree sequence $= (3, 3, 2, 2, 2, 2, 1, 1, 1, 1)$" + "\n"
               r"graph6: Ip\_I?D??G")
    )

    # n = 11: spider S^(2)_5 (hub degree 5)
    G11 = spider(5)
    pos11 = spider_layout(5)
    draw_tree(
        axes[1], G11, pos11,
        title=(r"$n = 11$, degree sequence $= (5, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1)$" + "\n"
               r"graph6: JqD?I?@O??\_")
    )

    # n = 12: DS(4, 3)
    G12 = double_spider(4, 3)
    pos12 = double_spider_layout(4, 3)
    draw_tree(
        axes[2], G12, pos12,
        title=(r"$n = 12$, degree sequence $= (4, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1)$" + "\n"
               r"graph6: Kp\_I?D??I??@")
    )

    # CRITICAL: the entire title is now plain text (NOT mathtext) with
    # explicit ASCII hyphen-minus characters. mathtext math-mode minus
    # signs (longer horizontal "U+2212"-style lines) were rendering thin
    # enough at the compiled-MATCH A5/10pt font size that a reviewer
    # mistook them for missing, parsing the title as W(T; 1, 1, 1).
    # Plain ASCII hyphen-minus rendered in bold serif is visually
    # unambiguous at every PDF zoom level.
    fig.suptitle(
        "Maximizers of W(T; -1, -1, 1) over trees of order n  (non-classical regime)",
        fontsize=14,
        weight="bold",
    )

    plt.tight_layout(rect=[0, 0, 1, 0.95])

    png_path = os.path.join(FIG_DIR, "fig10_caterpillar_maximizers.png")
    pdf_path = os.path.join(FIG_DIR, "fig10_caterpillar_maximizers.pdf")
    fig.savefig(png_path)
    fig.savefig(pdf_path)
    plt.close(fig)
    print(f"Wrote {png_path}")
    print(f"Wrote {pdf_path}")

    # Sanity-check the trees have the expected degree sequences
    for label, G, expected in [
        ("n=10 DS(3,3)", G10, sorted([3, 3, 2, 2, 2, 2, 1, 1, 1, 1], reverse=True)),
        ("n=11 spider 5", G11, sorted([5, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1], reverse=True)),
        ("n=12 DS(4,3)", G12, sorted([4, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1], reverse=True)),
    ]:
        actual = sorted([d for _, d in G.degree()], reverse=True)
        ok = actual == expected
        print(f"  {label:15s}  degrees={actual}  expected={expected}  {'OK' if ok else 'MISMATCH'}")


if __name__ == "__main__":
    main()
