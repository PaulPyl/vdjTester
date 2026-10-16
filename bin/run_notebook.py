#!/usr/bin/env python3
"""Execute a Jupyter notebook in this process, without a kernel, and save it.

Quarto normally executes Python through a Jupyter kernel. Under Rosetta (amd64
containers in the Mac's Lima VM) the kernel never answered its client, over TCP or
IPC, so `quarto render` timed out. This runs each code cell in one IPython
InteractiveShell instead and stores what the cell printed and displayed (text,
markdown, figures) as notebook outputs, in the order it happened, so a loop that
prints a heading and then shows a figure keeps them together. `quarto render` on
the result then only converts it to HTML; .ipynb files are not executed unless
asked.

Usage: run_notebook.py in.ipynb out.ipynb
"""

import contextlib
import sys

import matplotlib
import nbformat
from IPython.core.displaypub import DisplayPublisher
from IPython.core.interactiveshell import InteractiveShell


class Outputs:
    """Notebook outputs in arrival order; consecutive prints merge into one stream."""

    def __init__(self):
        self.items = []

    def write(self, text):
        if not text:
            return
        if self.items and self.items[-1].output_type == "stream":
            self.items[-1].text += text
        else:
            self.items.append(nbformat.v4.new_output("stream", name="stdout", text=text))

    def flush(self):
        pass

    def display(self, data, metadata):
        self.items.append(nbformat.v4.new_output("display_data", data=data, metadata=metadata or {}))


class OrderedPublisher(DisplayPublisher):
    """Sends display() calls, including inline figures, to the current Outputs."""

    target = None

    def publish(self, data, metadata=None, source=None, *, transient=None, update=False, **kwargs):
        if self.target is not None:
            self.target.display(data, metadata)


def main(src, dst):
    nb = nbformat.read(src, as_version=4)
    # The inline backend turns plt.show() into display(figure).
    # (%matplotlib inline needs a GUI hook the plain InteractiveShell lacks.)
    matplotlib.use("module://matplotlib_inline.backend_inline")
    shell = InteractiveShell.instance()
    publisher = OrderedPublisher(shell=shell)
    shell.display_pub = publisher
    for i, cell in enumerate(nb.cells):
        if cell.cell_type != "code":
            continue
        out = Outputs()
        publisher.target = out
        with contextlib.redirect_stdout(out):
            result = shell.run_cell(cell.source, store_history=False)
        publisher.target = None
        cell.outputs = out.items
        cell.execution_count = i + 1
        if not result.success:
            err = result.error_in_exec or result.error_before_exec
            sys.exit(f"cell {i + 1} failed: {type(err).__name__}: {err}\n--- cell source ---\n{cell.source}")
    nbformat.write(nb, dst)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
