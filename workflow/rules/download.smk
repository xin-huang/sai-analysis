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


rule download_lithuanian_genomes:
    input:
    output:
        map = "resources/data/Lithuanians/lit.map",
        ped = "resources/data/Lithuanians/lit.ped",
        download_flag = "resources/flags/.lithuanian_genomes.downloaded",
    shell:
        """
        wget -c https://data.mendeley.com/public-files/datasets/d2xt5hdm5j/files/04ec1cb6-f72f-4fc6-b443-db932cc59028/file_downloaded
        mv file_downloaded {output.map}
        wget -c https://data.mendeley.com/public-files/datasets/d2xt5hdm5j/files/ac649f46-2ef1-4650-9c7f-bf279de07d3f/file_downloaded
        mv file_downloaded {output.ped}
        touch {output.download_flag}
        """


rule download_nea_genome:
    input:
    output:
        vcf = "resources/data/NeaAltai/chr{i}_mq25_mapab100.vcf.gz",
        index = "resources/data/NeaAltai/chr{i}_mq25_mapab100.vcf.gz.tbi",
    shell:
        """
        wget -c http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Altai/chr{wildcards.i}_mq25_mapab100.vcf.gz
        mv chr{wildcards.i}_mq25_mapab100.vcf.gz {output.vcf}
        wget -c http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Altai/chr{wildcards.i}_mq25_mapab100.vcf.gz.tbi
        mv chr{wildcards.i}_mq25_mapab100.vcf.gz.tbi {output.index}
        """


rule download_den_genome:
    input:
    output:
        vcf = "resources/data/Denisova/chr{i}_mq25_mapab100.vcf.gz",
        index = "resources/data/Denisova/chr{i}_mq25_mapab100.vcf.gz.tbi",
    shell:
        """
        wget -c http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Denisova/chr{wildcards.i}_mq25_mapab100.vcf.gz
        mv chr{wildcards.i}_mq25_mapab100.vcf.gz {output.vcf}
        wget -c http://cdna.eva.mpg.de/neandertal/Vindija/VCF/Denisova/chr{wildcards.i}_mq25_mapab100.vcf.gz.tbi
        mv chr{wildcards.i}_mq25_mapab100.vcf.gz.tbi {output.index}
        """


rule download_1KG_genomes:
    input:
    output:
        vcf = "resources/data/1KG/ALL.chr{i}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
        index = "resources/data/1KG/ALL.chr{i}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi",
    shell:
        """
        wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{wildcards.i}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
        mv ALL.chr{wildcards.i}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz {output.vcf}
        wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{wildcards.i}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi
        mv ALL.chr{wildcards.i}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi {output.index}
        """


rule download_1KG_info:
    input:
    output:
        ref = "resources/data/1KG/human_g1k_v37.fasta",
        samples = "resources/data/1KG/integrated_call_samples_v3.20130502.ALL.panel",
    shell:
        """
        wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
        mv integrated_call_samples_v3.20130502.ALL.panel {output.samples}
        wget -c https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
        gzip -d human_g1k_v37.fasta.gz || true
        mv human_g1k_v37.fasta {output.ref}
        """


rule download_beagle:
    input:
    output:
        download_flag = "resources/flags/.beagle.downloaded",
    params:
        dir = "resources/tools/beagle/",
    shell:
        """
        mkdir -p {params.dir}
        cd {params.dir}
        wget -c https://faculty.washington.edu/browning/beagle/beagle.06Aug24.a91.jar
        cd ../../../
        touch {output.download_flag}
        """
