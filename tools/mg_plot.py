# -*- coding: utf-8 -*-
"""
Created on Sat Nov  6 15:31:33 2021

@author: mg
"""
# from collections.abc import Iterable

import matplotlib as mpl
import matplotlib.pyplot as plt

# mpl.rcParams["font.family"] = "Avenir LT Std"
plt.rcParams["font.size"] = 12
plt.rcParams["axes.linewidth"] = 1.5


def new_fig(*args, figsize=(8, 6), **kw_args):
    return plt.subplots(*args, figsize=figsize, **kw_args)


def set_styling(ax_or_axes, x_loc=None, y_loc=None, grid=False):
    if not isinstance(ax_or_axes, mpl.axes.Axes):
        for ax in ax_or_axes:
            set_styling(ax)
        return

    ax = ax_or_axes

    ax.xaxis.set_tick_params(which="major", size=10, width=1.5, direction="in", top="on")
    ax.xaxis.set_tick_params(which="minor", size=7, width=1.5, direction="in", top="on")
    ax.yaxis.set_tick_params(which="major", size=10, width=1.5, direction="in", right="on")
    ax.yaxis.set_tick_params(which="minor", size=7, width=1.5, direction="in", right="on")

    if x_loc is not None:
        if type(x_loc) is int or type(x_loc) is float:
            ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(x_loc))
        else:
            ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(x_loc[0]))
            ax.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(x_loc[1]))

    if y_loc is not None:
        if type(y_loc) is int or type(y_loc) is float:
            ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(y_loc))
        else:
            ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(y_loc[0]))
            ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(y_loc[1]))

    if grid:
        ax.grid(linestyle="dotted", linewidth=1.5)
