# Copyright 2025 Xin Huang
#
# GNU General Public License v3.0
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, please see
#
#    https://www.gnu.org/licenses/gpl-3.0.en.html


ruleorder: plot_1KG_outliers > plot_1KG_scores


rule plot_1KG_scores:
    input:
        u_scores = "results/sai/2src/1KG/nea_den/w_{w}_y_{y}_z_{z}/{allele_type}/1KG.nea_den.{approach}.w_{w}_y_{y}_z_{z}.U50.{allele_type}.scores.tsv",
        q_scores = "results/sai/2src/1KG/nea_den/w_{w}_y_{y}_z_{z}/{allele_type}/1KG.nea_den.{approach}.w_{w}_y_{y}_z_{z}.Q95.{allele_type}.scores.tsv",
    output:
        score_plot = "results/plots/2src/1KG/nea_den/w_{w}_y_{y}_z_{z}/{allele_type}/1KG.nea_den.{approach}.w_{w}_y_{y}_z_{z}.{allele_type}.scores.png",
    shell:
        """
        sai plot --u-file {input.u_scores} --q-file {input.q_scores} --title "" --xlabel Q95 --ylabel U50 --output {output.score_plot}
        """


rule plot_1KG_outliers:
    input:
        u_outliers = "results/sai/2src/1KG/nea_den/w_{w}_y_{y}_z_{z}/{allele_type}/1KG.nea_den.{approach}.w_{w}_y_{y}_z_{z}.U50.{allele_type}.outliers.tsv",
        q_outliers = "results/sai/2src/1KG/nea_den/w_{w}_y_{y}_z_{z}/{allele_type}/1KG.nea_den.{approach}.w_{w}_y_{y}_z_{z}.Q95.{allele_type}.outliers.tsv",
    output:
        outlier_plot = "results/plots/2src/1KG/nea_den/w_{w}_y_{y}_z_{z}/{allele_type}/1KG.nea_den.{approach}.w_{w}_y_{y}_z_{z}.{allele_type}.outliers.png",
        outlier_overlap = "results/plots/2src/1KG/nea_den/w_{w}_y_{y}_z_{z}/{allele_type}/1KG.nea_den.{approach}.w_{w}_y_{y}_z_{z}.{allele_type}.outliers.overlap.tsv",
    shell:
        """
        sai plot --u-file {input.u_outliers} --q-file {input.q_outliers} --title "" --xlabel Q95 --ylabel U50 --output {output.outlier_plot}
        """


rule plot_lit:
    input:
        u_scores = "results/sai/1src/lit/nea/w_{w}_y_{y}/{allele_type}/lit.nea.w_{w}_y_{y}.U90.{allele_type}.scores.tsv", 
        q_scores = "results/sai/1src/lit/nea/w_{w}_y_{y}/{allele_type}/lit.nea.w_{w}_y_{y}.Q95.{allele_type}.scores.tsv",
        u_outliers = "results/sai/1src/lit/nea/w_{w}_y_{y}/{allele_type}/lit.nea.w_{w}_y_{y}.U90.{allele_type}.outliers.tsv", 
        q_outliers = "results/sai/1src/lit/nea/w_{w}_y_{y}/{allele_type}/lit.nea.w_{w}_y_{y}.Q95.{allele_type}.outliers.tsv",
    output:
        score_plot = "results/plots/1src/lit/nea/w_{w}_y_{y}/{allele_type}/lit.nea.w_{w}_y_{y}.{allele_type}.scores.png",
        outlier_plot = "results/plots/1src/lit/nea/w_{w}_y_{y}/{allele_type}/lit.nea.w_{w}_y_{y}.{allele_type}.outliers.png",
        outlier_overlap = "results/plots/1src/lit/nea/w_{w}_y_{y}/{allele_type}/lit.nea.w_{w}_y_{y}.{allele_type}.outliers.overlap.tsv",
    shell:
        """
        sai plot --u-file {input.u_scores} --q-file {input.q_scores} --title "" --xlabel Q95 --ylabel U90 --output {output.score_plot}
        sai plot --u-file {input.u_outliers} --q-file {input.q_outliers} --title "" --xlabel Q95 --ylabel U90 --output {output.outlier_plot}
        """


rule plot_pan:
    input:
        u_scores = "results/sai/1src/pan/PPA/w_{w}_y_{y}/{allele_type}/pan.PPA.w_{w}_y_{y}.U90.{allele_type}.scores.tsv", 
        q_scores = "results/sai/1src/pan/PPA/w_{w}_y_{y}/{allele_type}/pan.PPA.w_{w}_y_{y}.Q95.{allele_type}.scores.tsv",
        u_outliers = "results/sai/1src/pan/PPA/w_{w}_y_{y}/{allele_type}/pan.PPA.w_{w}_y_{y}.U90.{allele_type}.outliers.tsv", 
        q_outliers = "results/sai/1src/pan/PPA/w_{w}_y_{y}/{allele_type}/pan.PPA.w_{w}_y_{y}.Q95.{allele_type}.outliers.tsv",
    output:
        score_plot = "results/plots/1src/pan/PPA/w_{w}_y_{y}/{allele_type}/pan.PPA.w_{w}_y_{y}.{allele_type}.scores.png",
        outlier_plot = "results/plots/1src/pan/PPA/w_{w}_y_{y}/{allele_type}/pan.PPA.w_{w}_y_{y}.{allele_type}.outliers.png",
        outlier_overlap = "results/plots/1src/pan/PPA/w_{w}_y_{y}/{allele_type}/pan.PPA.w_{w}_y_{y}.{allele_type}.outliers.overlap.tsv",
    shell:
        """
        sai plot --u-file {input.u_scores} --q-file {input.q_scores} --title "" --xlabel Q95 --ylabel U90 --output {output.score_plot}
        sai plot --u-file {input.u_outliers} --q-file {input.q_outliers} --title "" --xlabel Q95 --ylabel U90 --output {output.outlier_plot}
        """
