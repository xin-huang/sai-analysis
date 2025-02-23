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


rule annotate_1KG_candidates:
    input:
        outliers = "results/plots/2src/1KG/nea_den/w_{w}_y_{y}_z_{z}/{allele_type}/1KG.nea_den.{approach}.w_{w}_y_{y}_z_{z}.{allele_type}.outliers.overlap.tsv",
        annotation = expand("results/annotated_data/1KG/1KG.nea_den.{approach}.chr{i}.hg19_multianno.txt", i=list(range(1,23)), allow_missing=True),
    output:
        annotated_candidates = "results/plots/2src/1KG/nea_den/w_{w}_y_{y}_z_{z}/{allele_type}/1KG.nea_den.{approach}.w_{w}_y_{y}_z_{z}.{allele_type}.outliers.overlap.annotated.tsv",
    resources:
        mem_gb = 32,
    script:
       "../scripts/get_annotated_candidates.py"


rule annotate_lit_candidates:
    input:
        outliers = "results/plots/1src/lit/nea/w_{w}_y_{y}/{allele_type}/lit.nea.w_{w}_y_{y}.{allele_type}.outliers.overlap.tsv",
        annotation = expand("results/annotated_data/lit/lit.nea.chr{i}.hg19_multianno.txt", i=list(range(1,23)), allow_missing=True),
    output:
        annotated_candidates = "results/plots/1src/lit/nea/w_{w}_y_{y}/{allele_type}/lit.nea.w_{w}_y_{y}.{allele_type}.outliers.overlap.annotated.tsv",
    resources:
        mem_gb = 32,
    script:
       "../scripts/get_annotated_candidates.py"


rule annotate_pan_candidates:
    input:
        outliers = "results/plots/1src/pan/PPA/w_{w}_y_{y}/{allele_type}/pan.PPA.w_{w}_y_{y}.{allele_type}.{file_type}.overlap.tsv",
        annotation = expand("results/annotated_data/pan/pan.chr{i}.biallelic.snps.hg38_multianno.txt", i=list(range(1,23)), allow_missing=True),
    output:
        annotated_candidates = "results/plots/1src/pan/PPA/w_{w}_y_{y}/{allele_type}/pan.PPA.w_{w}_y_{y}.{allele_type}.{file_type}.overlap.annotated.tsv",
    resources:
        mem_gb = 32,
    script:
       "../scripts/get_annotated_candidates.py" 
