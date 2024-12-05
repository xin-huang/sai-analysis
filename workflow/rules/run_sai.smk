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


#rule analyze_1KG_samples:
#    input:
#        vcf = rules.merge_1KG_Nea.output.vcf,
#        ref = rules.create_ref_tgt_samples.output.ref,
#        tgt = rules.create_ref_tgt_samples.output.tgt,
#        src = rules.create_src_samples.output.nea_sample,
#    output:
#        scores = "results/sai/1KG/1KG.chr{i}.scores.txt",
#    params:
#        chr_name = "{i}",
#        win_len = 40000,
#        win_step = 40000,
#        w = 0.01,
#        x = 0.5,
#        y = 1,
#        q = 0.95,
#    resources:
#        cpus = 1, mem_gb = 48,
#    shell:
#        """
#        sai score --vcf {input.vcf} --ref {input.ref} --tgt {input.tgt} --src {input.src} --phased --w {params.w} --x {params.x} --y {params.y} --q {params.q} --chr-name {params.chr_name} --output {output.scores} --win-len {params.win_len} --win-step {params.win_step} --workers {resources.cpus}
#        """


rule analyze_lit_nea_samples:
    input:
        vcf = rules.merge_lit_yri_nea.output.vcf,
        ref = rules.extract_lit_samples.output.yri_samples,
        tgt = rules.extract_lit_samples.output.lit_samples,
        src = rules.extract_nea_sample.output.nea_sample,
    output:
        scores = "results/sai/Lithuanians/nea/w_{w}_x_{x}_y_{y}/lit.nea.chr{i}.w_{w}_x_{x}_y_{y}.scores.txt",
        u_outliers = "results/sai/Lithuanians/nea/w_{w}_x_{x}_y_{y}/lit.nea.chr{i}.w_{w}_x_{x}_y_{y}_U_outliers.tsv",
        q_outliers = "results/sai/Lithuanians/nea/w_{w}_x_{x}_y_{y}/lit.nea.chr{i}.w_{w}_x_{x}_y_{y}_Q95_outliers.tsv",
    params:
        chr_name = "{i}",
        win_len = 50000,
        win_step = 10000,
        w = "{w}",
        x = "{x}",
        y = "{y}",
        q = 0.95,
        output_dir = "results/sai/Lithuanians/nea/w_{w}_x_{x}_y_{y}",
        output_prefix = "lit.nea.chr{i}.w_{w}_x_{x}_y_{y}",
        outlier_quantile = 0.99,
    resources:
        cpus = 1, mem_gb = 32,
    shell:
        """
        sai score --vcf {input.vcf} --ref {input.ref} --tgt {input.tgt} --src {input.src} --phased --w {params.w} --x {params.x} --y {params.y} --q {params.q} --chr-name {params.chr_name} --output {output.scores} --win-len {params.win_len} --win-step {params.win_step} --workers {resources.cpus}
        sai outlier --score {output.scores} --output-dir {params.output_dir} --output-prefix {params.output_prefix} --quantile {params.outlier_quantile}
        """


rule analyze_lit_den_samples:
    input:
        vcf = rules.merge_lit_yri_den.output.vcf,
        ref = rules.extract_lit_samples.output.yri_samples,
        tgt = rules.extract_lit_samples.output.lit_samples,
        src = rules.extract_den_sample.output.den_sample,
    output:
        scores = "results/sai/Lithuanians/den/w_{w}_x_{x}_y_{y}/lit.den.chr{i}.w_{w}_x_{x}_y_{y}.scores.txt",
        u_outliers = "results/sai/Lithuanians/den/w_{w}_x_{x}_y_{y}/lit.den.chr{i}.w_{w}_x_{x}_y_{y}_U_outliers.tsv",
        q_outliers = "results/sai/Lithuanians/den/w_{w}_x_{x}_y_{y}/lit.den.chr{i}.w_{w}_x_{x}_y_{y}_Q95_outliers.tsv",
    params:
        chr_name = "{i}",
        win_len = 50000,
        win_step = 10000,
        w = "{w}",
        x = "{x}",
        y = "{y}",
        q = 0.95,
        output_dir = "results/sai/Lithuanians/den/w_{w}_x_{x}_y_{y}",
        output_prefix = "lit.den.chr{i}.w_{w}_x_{x}_y_{y}",
        outlier_quantile = 0.99,
    resources:
        cpus = 1, mem_gb = 32,
    shell:
        """
        sai score --vcf {input.vcf} --ref {input.ref} --tgt {input.tgt} --src {input.src} --phased --w {params.w} --x {params.x} --y {params.y} --q {params.q} --chr-name {params.chr_name} --output {output.scores} --win-len {params.win_len} --win-step {params.win_step} --workers {resources.cpus}
        sai outlier --score {output.scores} --output-dir {params.output_dir} --output-prefix {params.output_prefix} --quantile {params.outlier_quantile}
        """


rule analyze_lit_nea_den_samples:
    input:
        vcf = rules.merge_lit_yri_nea_den.output.vcf,
        ref = rules.extract_lit_samples.output.yri_samples,
        tgt = rules.extract_lit_samples.output.lit_samples,
        src = rules.extract_nea_den_samples.output.nea_den_samples,
    output:
        scores = "results/sai/Lithuanians/nea_den/w_{w}_x_{x}_y_{y}_z_{z}/lit.nea.den.chr{i}.w_{w}_x_{x}_y_{y}_z_{z}.scores.txt",
        u_outliers = "results/sai/Lithuanians/nea_den/w_{w}_x_{x}_y_{y}_z_{z}/lit.nea.den.chr{i}.w_{w}_x_{x}_y_{y}_z_{z}_U_outliers.tsv",
        q_outliers = "results/sai/Lithuanians/nea_den/w_{w}_x_{x}_y_{y}_z_{z}/lit.nea.den.chr{i}.w_{w}_x_{x}_y_{y}_z_{z}_Q95_outliers.tsv",
    params:
        chr_name = "{i}",
        win_len = 50000,
        win_step = 10000,
        w = "{w}",
        x = "{x}",
        y = "{y}",
        z = "{z}",
        q = 0.95,
        output_dir = "results/sai/Lithuanians/nea_den/w_{w}_x_{x}_y_{y}_z_{z}",
        output_prefix = "lit.nea.den.chr{i}.w_{w}_x_{x}_y_{y}_z_{z}",
        outlier_quantile = 0.99,
    resources:
        cpus = 1, mem_gb = 32,
    shell:
        """
        sai score --vcf {input.vcf} --ref {input.ref} --tgt {input.tgt} --src {input.src} --phased --w {params.w} --x {params.x} --y {params.y} --q {params.q} --chr-name {params.chr_name} --output {output.scores} --win-len {params.win_len} --win-step {params.win_step} --workers {resources.cpus} --num-src 2
        sai outlier --score {output.scores} --output-dir {params.output_dir} --output-prefix {params.output_prefix} --quantile {params.outlier_quantile}
        """
