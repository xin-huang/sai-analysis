# Copyright 2024 Xin Huang
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


rule analyze_1KG_2src_samples:
    input:
        vcf = rules.merge_1KG_nea_den.output.vcf,
        ref = rules.create_ref_tgt_samples.output.ref,
        tgt = rules.create_ref_tgt_samples.output.tgt,
        src = rules.extract_nea_den_samples.output.nea_den_samples,
    output:
        scores = "results/sai/2src/1KG/nea_den/w_{w}_x_{x}_y_{y}_z_{z}/1KG.nea_den.chr{i}.w_{w}_x_{x}_y_{y}_z_{z}.scores.tsv",
    params:
        chr_name = "{i}",
        win_len = 40000,
        win_step = 40000,
        w = "{w}",
        x = "{x}",
        y = "{y} {z}",
        q = 0.95,
    resources:
        cpus = 1, mem_gb = 48,
    shell:
        """
        sai score --vcf {input.vcf} --ref {input.ref} --tgt {input.tgt} --src {input.src} --phased --w {params.w} --x {params.x} --y {params.y} --q {params.q} --chr-name {params.chr_name} --output {output.scores} --win-len {params.win_len} --win-step {params.win_step} --workers {resources.cpus} --num-src 2
        """


rule analyze_lit_2src_samples:
    input:
        vcf = rules.merge_lit_yri_nea_den.output.vcf,
        ref = rules.extract_lit_samples.output.yri_samples,
        tgt = rules.extract_lit_samples.output.lit_samples,
        src = rules.extract_nea_den_samples.output.nea_den_samples,
    output:
        scores = "results/sai/2src/lit/nea_den/w_{w}_x_{x}_y_{y}_z_{z}/lit.nea_den.chr{i}.w_{w}_x_{x}_y_{y}_z_{z}.scores.tsv",
    params:
        chr_name = "{i}",
        win_len = 50000,
        win_step = 10000,
        w = "{w}",
        x = "{x}",
        y = "{y} {z}",
        q = 0.95,
    resources:
        cpus = 1, mem_gb = 32,
    shell:
        """
        sai score --vcf {input.vcf} --ref {input.ref} --tgt {input.tgt} --src {input.src} --phased --w {params.w} --x {params.x} --y {params.y} --q {params.q} --chr-name {params.chr_name} --output {output.scores} --win-len {params.win_len} --win-step {params.win_step} --workers {resources.cpus} --num-src 2
        """


rule get_2src_samples_outliers:
    input:
        scores = expand("results/sai/2src/{dataset}/nea_den/w_{w}_x_{x}_y_{y}_z_{z}/{dataset}.nea_den.chr{i}.w_{w}_x_{x}_y_{y}_z_{z}.scores.tsv", i=list(range(1,23)), allow_missing=True),
    output:
        all_scores = "results/sai/2src/{dataset}/nea_den/w_{w}_x_{x}_y_{y}_z_{z}/{dataset}.nea_den.w_{w}_x_{x}_y_{y}_z_{z}.scores.tsv",
        u_outliers = "results/sai/2src/{dataset}/nea_den/w_{w}_x_{x}_y_{y}_z_{z}/{dataset}.nea_den.w_{w}_x_{x}_y_{y}_z_{z}_U_outliers.tsv",
        q_outliers = "results/sai/2src/{dataset}/nea_den/w_{w}_x_{x}_y_{y}_z_{z}/{dataset}.nea_den.w_{w}_x_{x}_y_{y}_z_{z}_Q95_outliers.tsv",
    params:
        outlier_quantile = 0.99,
    shell:
        """
        set +o pipefail
        cat {input.scores} | head -1 > {output.all_scores}
        cat {input.scores} | grep -v "Chrom" >> {output.all_scores}
        sai outlier --score {output.all_scores} --output {output.u_outliers} --stat-type U --quantile {params.outlier_quantile}
        sai outlier --score {output.all_scores} --output {output.q_outliers} --stat-type Q --quantile {params.outlier_quantile}
        """