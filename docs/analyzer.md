---
template: main.html
title: gmx_MMPBSA_ana
---

# `gmx_MMPBSA_ana`: The analyzer tool

## Overview
`gmx_MMPBSA_ana` is a simple but powerful analysis tool. It is mainly focused on providing a fast, easy and efficient 
access to different graphics according to the analysis (**Figure 1**). The tool has been optimized to work with many 
charts, even more than you can review ([Check this section to learn about the potential of 
`gmx_MMPBSA_ana`](#gmx_mmpbsa_ana-under-pressure)).

<figure markdown="1">
[![overview][2]][2]
  <figcaption markdown="1" style="margin-top:0;">
**Figure 1.** `gmx_MMPBSA_ana` graphical overview
  </figcaption>
</figure>

[2]: assets/images/gmx_mmpbsa_ana_overview.png

**gmx_MMPBSA_ana components**

1. Start dialog
2. Plot area
3. Data control area
4. Menus

## Types of charts

### Line plot
Represents the evolution of a component during the simulation time. In some cases, it is the value of the calculated 
parameter, eg: `TOTAL DELTA`, `VDWAALS`, etc., and in others, it is the sum of the elements it contains, eg: `TDC` as 
the sum of all the per-residue energy contributions.

!!! note
    The line plot can additionally contain the moving average (solid red line) or other elements such as other 
    indicators (dash red and green lines).

!!! warning 
    Note that `TDC`, `SDC`, and `BDC` are calculated from the immediate items contained within themselves. This means that the 
    results obtained for the `TDC` in a per-residue calculation will probably be different from the `TDC` of a per-wise 
    calculation with the same residues. This is because the per-residue energy contribution is calculated taking 
    into account the environment that surrounds each residue, so the sum of their contributions will be exactly the 
    maximum contribution per selected residues. On the other hand, in the per-wise calculation, the `TDC` will be 
    equal to the sum of each residue contribution, the difference being that the contribution of each residue is 
    obtained from the contributions of the selected pairs and not with all the environment that surrounds it. Therefore, both 
    calculations should be used for different purposes and care should be taken with the interpretation of the results.

<figure markdown="1">
[![lineplot][3]][3]
[![lineplot2][4]][4]
  <figcaption markdown="1" style="margin-top:0;">
  **Figure 2**. Line plot examples. **Up:** ΔH representation, **Down:** Interaction Entropy representation
  </figcaption>
</figure>

[3]: assets/images/line_plot.png
[4]: assets/images/line_plot_ie.png

### Bar plot
Represents the total contribution of a component for the simulation. In some cases, it is the Average of the calculated 
parameter, eg: `TOTAL DELTA`, `VDWAALS`, RESIDUES in per-residue calculation, etc., and in others, it is the sum of the 
elements it contains, eg: NMODE and QH Entropy, ΔG Binding, etc. The bars that represent averages also have a solid 
line that represents the standard deviation, while the bars representing a specific value do not.

<figure markdown="1">
[![barplot][5]][5]
[![barplot2][6]][6]
  <figcaption markdown="1" style="margin-top:0;">
  **Figure 3**. Bar plot examples. **Up:** Per-residue contribution, **Down:** ΔG Binding
  </figcaption>
</figure>

[5]: assets/images/bar_plot.png
[6]: assets/images/bar_plot2.png

### Heatmap plot
Represents the evolution of several components at the same time during the simulation. It is used to represent the 
contribution of all residues in per-residue calculation or all pairs for each residue in the per-wise calculation. 
It could also represent the relationship between components.This is explicit for per-wise calculations and 
represents the contribution of each residue with its respective pairs.

!!! tip
    The relational heatmap graph is the best representation in the per-wise analysis.

<figure markdown="1">
[![heatmapplot][7]][7]
[![heatmapplot2][8]][8]
  <figcaption markdown="1" style="margin-top:0;">
  **Figure 4**. Heatmap plot examples. **Up:** Per-residue contribution per-frame, **Down:** Inter Residue-pair 
contribution
  </figcaption>
</figure>

[7]: assets/images/heatmap_plot.png
[8]: assets/images/heatmap_plot2.png

---------------------------------------

## Representations

!!! note
    The following videos provide a visual medium to facilitate the use of `gmx_MMPBSA_ana`. If you find our tool 
    useful and want to use it for a visual project, such as tutorials, examples, etc., we will gladly include it 
    in the documentation as well as in the list of acknowledgments and collaborators.

### Binding Free Energy calculation (GB + IE) 

<div class="embed-container">
    <iframe src="https://www.youtube.com/embed/k1aLlBhnkxo" frameborder="0" allowfullscreen></iframe>
</div>

### Correlation analysis

<div class="embed-container">
    <iframe src="https://www.youtube.com/embed/0xiphzA1O0w" frameborder="0" allowfullscreen></iframe>
</div>

### Alanine scanning

<div class="embed-container">
    <iframe src="https://www.youtube.com/embed/_r13tcmY038" frameborder="0" allowfullscreen></iframe>
</div>

### Per-residue decomposition

<div class="embed-container">
    <iframe src="https://www.youtube.com/embed/Ww7juWeWQQ8" frameborder="0" allowfullscreen></iframe>
</div>

### Functionalities

<div class="embed-container">
    <iframe src="https://www.youtube.com/embed/PgDnG8UgRWw" frameborder="0" allowfullscreen></iframe>
</div>

## `gmx_MMPBSA_ana` under pressure

!!! warning
    Around 1.8 million graphs were loaded in gmx_MMPBSA_ana, and it is showed only for educational purposes.

<div class="embed-container">
    <iframe src="https://www.youtube.com/embed/xRIi3LtB7wU" frameborder="0" allowfullscreen></iframe>
</div>

In this experiment, we replicated our examples folder 9 times (differs from the one available on GitHub) giving us a 
total of 99 systems. As you can see in the video, gmx_MMPBSA_ana manages to deal well with the incredible amount of ~1.6 
million items. Every item contains between 1 and 3 graphics, for a total of ~1.8 million graphics loaded. This feat is 
accomplished in ~11 minutes. Most of this time is consumed processing the data associated with each graph. At the 
moment, gmx_MMPBSA_ana processes this data serially, since parallelizing this process would be a bit difficult. In 
any case, for the usual processes, this will take a maximum of 25-30 seconds, depending on your hardware. Each item has 
associated the data of each of its graphs, which is stored in memory. In this experiment, RAM consumption reached up 
to 14 GB.

!!! danger
    Be aware that if you run out of available RAM, your OS could crash, freeze, or slow down.

Does this mean that you will not be able to load 100 systems in `gmx_MMPBSA_ana`?
Not at all. Consumption depends on the type of calculation you have made, and the data you want to analyze. In this
experiment, there were several systems that contain the `decomp` data. We also selected to show the complex, receptor, and 
ligand data, which usually can be skipped. Only a system calculated with per-wise selecting about 40 amino acids 
generates about 11 thousand items.

